"""Unit tests for aaftf/pipeline.py: every step gets its own CLI defaults plus the pipeline's options."""

import argparse as ap
from pathlib import Path
from unittest.mock import patch

import pytest

from aaftf import pipeline
from aaftf._menu import register_subcommands

pytestmark = pytest.mark.unit

STEPS = ["trim", "filter", "assemble", "vecscreen", "taxpurge", "covpurge", "rmdup", "sort", "assess"]
_MODULE = {"filter": "aaftf_filter", "sort": "aaftf_sort"}
# the kwarg that names each step's output file (trim and filter derive theirs from basename)
_OUTPUT = {"assemble": "out", "polish": "outfile", "vecscreen": "outfile", "fcs_screen": "outfile", "fcs_gx_purge": "outfile", "taxpurge": "outfile", "covpurge": "outfile", "rmdup": "out", "sort": "out"}


def _cli_defaults(step):
    """Return ``{dest: default}`` for ``AAFTF <step>``."""
    parser = register_subcommands(ap.ArgumentParser()).choices[step]
    return {a.dest: a.default for a in parser._actions if a.dest != "help"}


def _run_steps(func, tmp_path, **options):
    """Run a pipeline function with every step's run() replaced; return [(step, kwargs)] in call order."""
    base = str(tmp_path / "sample")
    calls = []

    def _fake(step):
        def _run(**kwargs):
            calls.append((step, kwargs))
            # trim writes _1P/_2P for paired reads and _1U for single-end reads, as aaftf.trim does
            trimmed = f"{base}_1P.fastq.gz" if kwargs.get("read2") else f"{base}_1U.fastq.gz"
            outputs = {"trim": [trimmed], "filter": [f"{base}_filtered_1.fastq.gz"]}.get(step, [kwargs.get(_OUTPUT.get(step, ""))])
            for out in filter(None, outputs):
                Path(out).write_text("data")

        return _run

    patches = [patch(f"aaftf.pipeline.{_MODULE.get(step, step)}.run", side_effect=_fake(step)) for step in STEPS + ["polish", "fcs_screen", "fcs_gx_purge"]]
    for p in patches:
        p.start()
    try:
        func(basename=base, **options)
    finally:
        for p in patches:
            p.stop()
    return calls


def _run_pipeline(tmp_path, read2=True, **options):
    """Run pipeline.run_short() with every step replaced; return {step: kwargs it was called with}."""
    return dict(_run_steps(pipeline.run_short, tmp_path, read1="R1.fq.gz", read2="R2.fq.gz" if read2 else None, phylum=["Ascomycota"], **options))


class TestSteps:
    def test_runs_every_step_and_no_optional_ones(self, tmp_path):
        assert list(_run_pipeline(tmp_path)) == STEPS

    @pytest.mark.parametrize("step", STEPS)
    def test_step_gets_every_cli_option(self, tmp_path, step):
        assert set(_run_pipeline(tmp_path)[step]) == set(_cli_defaults(step)) | {"pipe"}

    @pytest.mark.parametrize("step", STEPS)
    def test_step_hints_suppressed(self, tmp_path, step):
        assert _run_pipeline(tmp_path)[step]["pipe"] is True

    def test_existing_output_skips_step(self, tmp_path):
        Path(f"{tmp_path / 'sample'}_1P.fastq.gz").write_text("done")
        assert "trim" not in _run_pipeline(tmp_path)

    def test_missing_output_raises(self, tmp_path):
        with patch("aaftf.pipeline.trim.run"):
            with pytest.raises(RuntimeError, match="AAFTF trim failed"):
                pipeline.run_short(read1="R1.fq.gz", basename=str(tmp_path / "sample"), phylum=["Ascomycota"])


class TestDefaults:
    """Options the pipeline does not expose keep each step's own default."""

    @pytest.mark.parametrize(
        "step, option",
        [
            ("trim", "method"),
            ("trim", "avgqual"),
            ("filter", "aligner"),
            ("filter", "screen_local"),
            ("assemble", "isolate"),
            ("assemble", "careful"),
            ("assemble", "merged"),
            ("vecscreen", "percent_id"),
            ("vecscreen", "stringency"),
            ("taxpurge", "kmer"),
            ("taxpurge", "sourdb"),
            ("taxpurge", "min_hashes"),
            ("rmdup", "percent_id"),
            ("rmdup", "percent_cov"),
            ("sort", "name"),
            ("assess", "telomere_monomer"),
            ("assess", "report"),
        ],
    )
    def test_step_default_kept(self, tmp_path, step, option):
        assert _run_pipeline(tmp_path)[step][option] == _cli_defaults(step)[option]

    @pytest.mark.parametrize("step", ["trim", "filter", "assemble"])
    def test_memory_default_kept_without_pipeline_memory(self, tmp_path, step):
        assert _run_pipeline(tmp_path)[step]["memory"] == _cli_defaults(step)["memory"]


class TestPipelineOptions:
    """Options given to the pipeline override the step defaults."""

    def test_memory_passed_to_each_step(self, tmp_path):
        calls = _run_pipeline(tmp_path, memory=12)
        assert [calls[s]["memory"] for s in ("trim", "filter", "assemble")] == [12, 12, 12]

    @pytest.mark.parametrize("step", ["trim", "filter", "assemble", "vecscreen", "covpurge", "rmdup"])
    def test_cpus_passed(self, tmp_path, step):
        assert _run_pipeline(tmp_path, cpus=6)[step]["cpus"] == 6

    @pytest.mark.parametrize("step", ["filter", "vecscreen", "taxpurge", "covpurge", "rmdup"])
    def test_workdir_passed(self, tmp_path, step):
        assert _run_pipeline(tmp_path, workdir="wd")[step]["workdir"] == "wd"

    @pytest.mark.parametrize("method", ["spades", "megahit", "unicycler"])
    def test_assemble_gets_own_workdir_subfolder(self, tmp_path, method):
        """filter creates the pipeline workdir first; SPAdes would restart from it and MEGAHIT refuse it."""
        assert _run_pipeline(tmp_path, workdir="wd", method=method)["assemble"]["workdir"] == str(Path("wd", f"assemble_{method}"))

    def test_assemble_workdir_default_without_pipeline_workdir(self, tmp_path):
        assert _run_pipeline(tmp_path)["assemble"]["workdir"] == _cli_defaults("assemble")["workdir"]

    @pytest.mark.parametrize("step", STEPS)
    def test_debug_passed(self, tmp_path, step):
        assert _run_pipeline(tmp_path, debug=True)[step]["debug"] is True

    def test_step_specific_options(self, tmp_path):
        calls = _run_pipeline(tmp_path, minlen=50, mincontiglen=1000, method="megahit", mincovpct=7, sourdb=["gbk-fungi"], screen_accessions=["NC_1"])
        assert calls["trim"]["minlen"] == 50
        assert calls["rmdup"]["minlen"] == calls["sort"]["minlen"] == 1000
        assert calls["assemble"]["method"] == "megahit"
        assert (calls["covpurge"]["mincovpct"], calls["taxpurge"]["sourdb"]) == (7, ["gbk-fungi"])
        assert calls["taxpurge"]["phylum"] == ["Ascomycota"]
        assert calls["filter"]["screen_accessions"] == ["NC_1"]

    def test_pipeline_method_is_only_the_assembler(self, tmp_path):
        calls = _run_pipeline(tmp_path, method="megahit")
        assert calls["trim"]["method"] == _cli_defaults("trim")["method"]


class TestFileChaining:
    def test_each_step_reads_the_previous_output(self, tmp_path):
        calls = _run_pipeline(tmp_path)
        base = str(tmp_path / "sample")
        assert (calls["filter"]["read1"], calls["filter"]["read2"]) == (f"{base}_1P.fastq.gz", f"{base}_2P.fastq.gz")
        assert calls["assemble"]["read1"] == f"{base}_filtered_1.fastq.gz"
        assert calls["vecscreen"]["infile"] == f"{base}.spades.fasta"
        assert calls["taxpurge"]["input"] == f"{base}.vecscreen.fasta"
        assert calls["covpurge"]["input"] == f"{base}.taxpurge.fasta"
        assert (calls["covpurge"]["read1"], calls["covpurge"]["read2"]) == (f"{base}_filtered_1.fastq.gz", f"{base}_filtered_2.fastq.gz")
        assert calls["rmdup"]["input"] == f"{base}.covpurge.fasta"
        assert calls["sort"]["input"] == f"{base}.rmdup.fasta"  # polish is optional and not run
        assert calls["assess"]["input"] == f"{base}.final.fasta"

    def test_single_end_reads_have_no_read2(self, tmp_path):
        calls = _run_pipeline(tmp_path, read2=False)
        assert all(calls[s]["read2"] is None for s in ("trim", "filter", "assemble", "covpurge"))
        base = str(tmp_path / "sample")
        assert calls["filter"]["read1"] == f"{base}_1U.fastq.gz"  # trim's single-end output
        assert calls["assemble"]["read1"] == f"{base}_filtered_1.fastq.gz"

    def test_no_phylum_skips_taxpurge(self, tmp_path):
        base = str(tmp_path / "sample")
        calls = dict(_run_steps(pipeline.run_short, tmp_path, read1="R1.fq.gz", read2="R2.fq.gz"))
        assert list(calls) == [s for s in STEPS if s != "taxpurge"]
        assert calls["covpurge"]["input"] == f"{base}.vecscreen.fasta"


LONG_STEPS = ["assemble", "vecscreen", "rmdup", "polish", "sort", "assess"]
HYBRID_STEPS = ["trim", "filter", "assemble", "vecscreen", "taxpurge", "covpurge", "rmdup", "polish", "polish", "sort", "assess"]


def _run_long(tmp_path, **options):
    return _run_steps(pipeline.run_long, tmp_path, longreads="ont.fq.gz", **options)


def _run_hybrid(tmp_path, read2=True, **options):
    return _run_steps(pipeline.run_hybrid, tmp_path, read1="R1.fq.gz", read2="R2.fq.gz" if read2 else None, longreads="ont.fq.gz", phylum=["Ascomycota"], **options)


class TestLongPipeline:
    def test_step_order(self, tmp_path):
        assert [step for step, _ in _run_long(tmp_path)] == LONG_STEPS

    def test_every_step_gets_every_cli_option(self, tmp_path):
        for step, kwargs in _run_long(tmp_path):
            assert set(kwargs) == set(_cli_defaults(step)) | {"pipe"}, step

    def test_flye_then_racon_chain(self, tmp_path):
        calls = dict(_run_long(tmp_path, longread_type="pacbio-hifi", genome_size="40m", cpus=4))
        base = str(tmp_path / "sample")
        asm = calls["assemble"]
        assert (asm["method"], asm["longreads"], asm["read1"], asm["out"]) == ("flye", "ont.fq.gz", None, f"{base}.flye.fasta")
        assert (asm["longread_type"], asm["genome_size"], asm["cpus"]) == ("pacbio-hifi", "40m", 4)
        assert calls["vecscreen"]["infile"] == f"{base}.flye.fasta"
        assert calls["rmdup"]["input"] == f"{base}.vecscreen.fasta"  # no taxpurge without --phylum, no covpurge without Illumina reads
        assert (calls["polish"]["method"], calls["polish"]["infile"], calls["polish"]["longreads"]) == ("racon", f"{base}.rmdup.fasta", "ont.fq.gz")
        assert calls["sort"]["input"] == f"{base}.racon.fasta"

    def test_assemble_gets_own_workdir_subfolder(self, tmp_path):
        calls = dict(_run_long(tmp_path, workdir="wd"))
        assert (calls["assemble"]["workdir"], calls["polish"]["workdir"]) == (str(Path("wd", "assemble_flye")), "wd")

    def test_phylum_runs_taxpurge(self, tmp_path):
        calls = _run_long(tmp_path, phylum=["Ascomycota"], sourdb=["gbk-fungi"])
        base = str(tmp_path / "sample")
        assert [step for step, _ in calls] == ["assemble", "vecscreen", "taxpurge", "rmdup", "polish", "sort", "assess"]
        steps = dict(calls)
        purge = steps["taxpurge"]
        assert (purge["input"], purge["phylum"], purge["sourdb"]) == (f"{base}.vecscreen.fasta", ["Ascomycota"], ["gbk-fungi"])
        assert set(purge) == set(_cli_defaults("taxpurge")) | {"pipe"}
        assert steps["rmdup"]["input"] == f"{base}.taxpurge.fasta"

    def test_assemble_defaults_kept(self, tmp_path):
        asm = dict(_run_long(tmp_path))["assemble"]
        assert (asm["longread_type"], asm["genome_size"]) == (_cli_defaults("assemble")["longread_type"], None)


class TestHybridPipeline:
    def test_step_order(self, tmp_path):
        assert [step for step, _ in _run_hybrid(tmp_path)] == HYBRID_STEPS

    def test_every_step_gets_every_cli_option(self, tmp_path):
        for step, kwargs in _run_hybrid(tmp_path):
            assert set(kwargs) == set(_cli_defaults(step)) | {"pipe"}, step

    def test_flye_racon_polypolish_chain(self, tmp_path):
        calls = _run_hybrid(tmp_path)
        base = str(tmp_path / "sample")
        racon, short = [kwargs for step, kwargs in calls if step == "polish"]
        assert (racon["method"], racon["infile"]) == ("racon", f"{base}.rmdup.fasta")
        assert (short["method"], short["infile"]) == ("polypolish", f"{base}.racon.fasta")
        assert (short["read1"], short["read2"]) == (f"{base}_filtered_1.fastq.gz", f"{base}_filtered_2.fastq.gz")
        steps = dict(calls)
        assert steps["vecscreen"]["infile"] == f"{base}.flye.fasta"
        assert steps["sort"]["input"] == f"{base}.polypolish.fasta"
        assert steps["covpurge"]["read1"] == f"{base}_filtered_1.fastq.gz"

    def test_no_phylum_skips_taxpurge(self, tmp_path):
        calls = _run_steps(pipeline.run_hybrid, tmp_path, read1="R1.fq.gz", read2="R2.fq.gz", longreads="ont.fq.gz")
        assert [step for step, _ in calls] == [s for s in HYBRID_STEPS if s != "taxpurge"]

    def test_single_end_uses_pypolca(self, tmp_path):
        short = [kwargs for step, kwargs in _run_hybrid(tmp_path, read2=False) if step == "polish"][-1]
        assert (short["method"], short["read2"]) == ("pypolca", None)

    def test_unicycler_gets_both_read_types_and_skips_polish(self, tmp_path):
        calls = _run_hybrid(tmp_path, method="unicycler")
        base = str(tmp_path / "sample")
        assert "polish" not in [step for step, _ in calls]
        asm = dict(calls)["assemble"]
        assert (asm["method"], asm["read1"], asm["longreads"]) == ("unicycler", f"{base}_filtered_1.fastq.gz", "ont.fq.gz")
        assert dict(calls)["vecscreen"]["infile"] == f"{base}.unicycler.fasta"
        assert dict(calls)["sort"]["input"] == f"{base}.rmdup.fasta"

    @pytest.mark.parametrize("method", ["flye", "unicycler"])
    def test_assemble_gets_own_workdir_subfolder(self, tmp_path, method):
        steps = dict(_run_hybrid(tmp_path, workdir="wd", method=method))
        assert (steps["assemble"]["workdir"], steps["filter"]["workdir"]) == (str(Path("wd", f"assemble_{method}")), "wd")

    def test_unknown_method_raises(self, tmp_path):
        with pytest.raises(ValueError, match="flye or unicycler"):
            _run_hybrid(tmp_path, method="spades")


class TestPipelineParsers:
    """Each pipeline subcommand dispatches to its own run function with its own options."""

    @pytest.mark.parametrize(
        "argv, func, expected",
        [
            (["pipeline_short", "-1", "R1.fq", "-o", "s", "-p", "Ascomycota"], "run_short", {"read1": "R1.fq", "method": "spades"}),
            (["pipeline_long", "-lr", "ont.fq", "-o", "s", "--genome_size", "40m"], "run_long", {"longreads": "ont.fq", "genome_size": "40m", "longread_type": "nano-hq"}),
            (["pipeline_hybrid", "-1", "R1.fq", "-lr", "ont.fq", "-o", "s", "-p", "Ascomycota"], "run_hybrid", {"read1": "R1.fq", "longreads": "ont.fq", "method": "flye"}),
        ],
    )
    def test_dispatch(self, argv, func, expected):
        import sys

        from aaftf.main import main

        with patch.object(sys, "argv", ["AAFTF", *argv]), patch(f"aaftf.pipeline.{func}") as mock_run:
            # main() builds the parsers, binding the patched function, inside the patch
            assert main() == 0
        kwargs = mock_run.call_args.kwargs
        assert {k: kwargs[k] for k in expected} == expected

    def test_long_requires_longreads(self):
        import sys

        from aaftf.main import main

        with patch.object(sys, "argv", ["AAFTF", "pipeline_long", "-o", "s"]), pytest.raises(SystemExit):
            main()


class TestFcsSteps:
    """--fcs_screen / --fcs_gx_db add NCBI FCS-adaptor and FCS-GX right after vecscreen, in every pipeline."""

    def test_not_run_by_default(self, tmp_path):
        assert not {"fcs_screen", "fcs_gx_purge"} & set(_run_pipeline(tmp_path))

    def test_fcs_gx_replaces_taxpurge(self, tmp_path):
        """FCS-GX is the taxonomic screen instead of taxpurge: vecscreen → fcs_screen → fcs_gx_purge → covpurge."""
        base = str(tmp_path / "sample")
        calls = _run_steps(pipeline.run_short, tmp_path, read1="R1.fq.gz", read2="R2.fq.gz", fcs_screen_run=True, container_engine="docker", fcs_gx_db="/gxdb/all", taxid=5062)
        assert [step for step, _ in calls][3:7] == ["vecscreen", "fcs_screen", "fcs_gx_purge", "covpurge"]
        steps = dict(calls)
        assert "taxpurge" not in steps
        assert (steps["fcs_screen"]["infile"], steps["fcs_screen"]["container_engine"]) == (f"{base}.vecscreen.fasta", "docker")
        assert (steps["fcs_gx_purge"]["input"], steps["fcs_gx_purge"]["db"], steps["fcs_gx_purge"]["taxid"]) == (f"{base}.fcs_screen.fasta", "/gxdb/all", 5062)
        assert steps["covpurge"]["input"] == f"{base}.fcs_gx.fasta"
        for step in ("fcs_screen", "fcs_gx_purge"):
            assert set(steps[step]) == set(_cli_defaults(step)) | {"pipe"}, step

    def test_fcs_screen_with_taxpurge(self, tmp_path):
        """fcs_screen (adaptors, not taxonomy) combines with taxpurge."""
        base = str(tmp_path / "sample")
        steps = dict(_run_steps(pipeline.run_short, tmp_path, read1="R1.fq.gz", read2="R2.fq.gz", phylum=["Ascomycota"], fcs_screen_run=True))
        assert steps["taxpurge"]["input"] == f"{base}.fcs_screen.fasta" and "fcs_gx_purge" not in steps

    @pytest.mark.parametrize("func", [pipeline.run_short, pipeline.run_long, pipeline.run_hybrid])
    def test_phylum_and_fcs_gx_db_are_alternatives(self, tmp_path, func):
        with pytest.raises(ValueError, match="choose one taxonomic screen"):
            func(read1="R1.fq.gz", longreads="ont.fq.gz", basename=str(tmp_path / "s"), phylum=["Ascomycota"], fcs_gx_db="/gxdb/all", taxid=4890)

    def test_fcs_screen_alone_keeps_its_engine_default(self, tmp_path):
        steps = dict(_run_steps(pipeline.run_long, tmp_path, longreads="ont.fq.gz", fcs_screen_run=True))
        assert steps["fcs_screen"]["container_engine"] == _cli_defaults("fcs_screen")["container_engine"]
        assert steps["rmdup"]["input"].endswith(".fcs_screen.fasta")

    @pytest.mark.parametrize("func", [pipeline.run_short, pipeline.run_long, pipeline.run_hybrid])
    def test_fcs_gx_db_needs_taxid(self, tmp_path, func):
        with pytest.raises(ValueError, match="--taxid"):
            func(read1="R1.fq.gz", longreads="ont.fq.gz", basename=str(tmp_path / "s"), fcs_gx_db="/gxdb/all")

    def test_hybrid_runs_them_too(self, tmp_path):
        calls = _run_steps(pipeline.run_hybrid, tmp_path, read1="R1.fq.gz", read2="R2.fq.gz", longreads="ont.fq.gz", fcs_gx_db="/gxdb/all", taxid=4890)
        steps = [step for step, _ in calls]
        assert steps.index("vecscreen") < steps.index("fcs_gx_purge") < steps.index("covpurge")
