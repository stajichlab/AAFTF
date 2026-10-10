"""Unit tests for aaftf/mito.py helpers."""

from pathlib import Path
from unittest.mock import patch

import pytest

from aaftf.mito import run

pytestmark = pytest.mark.unit


class TestOrientToStart:
    """_orient_to_start rotates the circular genome so the start gene's alignment begins at position 0."""

    def _orient(self, tmp_path, hits, seq):
        from aaftf.mito import _orient_to_start
        from aaftf.utility import PafHit

        fasta_in, fasta_out = tmp_path / "in.fa", tmp_path / "out.fa"
        fasta_in.write_text(f">mt\n{seq}\n")
        paf = [PafHit("COB", 100, q_start, 100, strand, "mt", len(seq), t_start, t_end, 100, 100, 60) for q_start, strand, t_start, t_end in hits]
        with patch("aaftf.mito.paf_hits", return_value=iter(paf)):
            _orient_to_start(str(fasta_in), str(fasta_out))
        return "".join(fasta_out.read_text().splitlines()[1:])

    def test_forward_hit_rotates_to_target_start(self, tmp_path):
        seq = "AAAACCCCGGGGTTTT"
        assert self._orient(tmp_path, [(0, "+", 4, 8)], seq) == "CCCCGGGGTTTTAAAA"

    def test_reverse_hit_rotates_and_reverse_complements(self, tmp_path):
        seq = "AAAACCCCGGGGTTTA"  # not its own reverse complement once rotated
        # rotated to start at 8 -> GGGGTTTAAAAACCCC, then reverse complemented
        assert self._orient(tmp_path, [(0, "-", 4, 8)], seq) == "GGGGTTTTTAAACCCC"

    def test_no_hit_leaves_sequence_unrotated(self, tmp_path):
        seq = "AAAACCCCGGGGTTTT"
        assert self._orient(tmp_path, [], seq) == seq

    def test_negative_offset_leaves_sequence_unrotated(self, tmp_path):
        seq = "AAAACCCCGGGGTTTT"
        # query_start 10 pushes ref_start to 4 - 10 = -6 (would wrap via negative slicing)
        assert self._orient(tmp_path, [(10, "+", 4, 8)], seq) == seq

    def test_offset_past_end_leaves_sequence_unrotated(self, tmp_path):
        seq = "AAAACCCCGGGGTTTA"
        # minus strand: ref_start = 12 + 10 = 22 > len(seq); no rotation or reverse complement
        assert self._orient(tmp_path, [(10, "-", 8, 12)], seq) == seq

    def test_two_hits_leave_sequence_unrotated(self, tmp_path):
        seq = "AAAACCCCGGGGTTTT"
        assert self._orient(tmp_path, [(0, "+", 4, 8), (0, "+", 8, 12)], seq) == seq

    def _run(self, tmp_path, **kwargs):
        from aaftf.mito import run

        workdir = tmp_path / "wd"
        with patch("aaftf.mito.require_tools"), patch("aaftf.mito.estimate_read_length", return_value=150):
            with patch("aaftf.mito.subprocess.Popen") as popen:
                popen.return_value.communicate.return_value = (b"", b"")
                # no NOVOPlasty output is faked, so run() stops after writing its inputs
                with pytest.raises(RuntimeError, match="did not produce an assembly"):
                    run(read1=str(tmp_path / "R1.fq"), read2=str(tmp_path / "R2.fq"), out=str(tmp_path / "mt.fa"), workdir=str(workdir), subsample=0, **kwargs)
        return workdir

    def test_default_seed_copied_from_package_and_config_filled(self, tmp_path):
        workdir = self._run(tmp_path, memory=6)
        seed = workdir / "mito-seed.fasta"
        assert seed.read_text().startswith(">")
        config = (workdir / "novo-config.txt").read_text()
        assert "<" not in config.replace("<=", "")  # every <PLACEHOLDER> replaced
        assert str(seed.resolve()) in config and "150" in config and "6" in config

    def test_explicit_seed_not_copied(self, tmp_path):
        own_seed = tmp_path / "myseed.fa"
        own_seed.write_text(">s\nACGT\n")
        workdir = self._run(tmp_path, seed=str(own_seed))
        assert not (workdir / "mito-seed.fasta").exists()
        assert str(own_seed.resolve()) in (workdir / "novo-config.txt").read_text()


class TestNovoplastyOutputChoice:
    """run() picks NOVOPlasty's circular assembly first, then Contigs_1, then the uncircularised file."""

    def _run(self, tmp_path, outputs):
        def _fake_popen(cmd, cwd, **kwargs):
            for name, seq in outputs.items():
                Path(cwd, name).write_text(f">x\n{seq}\n")
            return type("P", (), {"communicate": lambda self: None})()

        orient = []
        out = tmp_path / "mt.fasta"
        with patch("aaftf.mito.require_tools"), patch("aaftf.mito.estimate_read_length", return_value=150):
            with patch("aaftf.mito.subprocess.Popen", side_effect=_fake_popen):
                with patch("aaftf.mito._orient_to_start", side_effect=lambda draft, out, **kw: orient.append(Path(draft).name)):
                    run(read1="R1.fq", read2="R2.fq", out=str(out), workdir=str(tmp_path / "wd"), subsample=0)
        return orient, out

    def test_circular_preferred(self, tmp_path):
        orient, _ = self._run(tmp_path, {"Contigs_1_x.fasta": "AAAA", "Circularized_assembly_1_x.fasta": "CCCC", "Uncircularized_assemblies_x.fasta": "GGGG"})
        assert orient == ["Circularized_assembly_1_x.fasta"]

    def test_contigs_preferred_over_uncircularized(self, tmp_path):
        orient, out = self._run(tmp_path, {"Uncircularized_assemblies_x.fasta": "GGGG", "Contigs_1_x.fasta": "AAAA"})
        assert orient == [] and "AAAA" in out.read_text()

    def test_no_assembly_raises(self, tmp_path):
        with pytest.raises(RuntimeError, match="NOVOplasty did not produce an assembly"):
            self._run(tmp_path, {})


class TestStartSequence:
    def test_bundled_cob_consensus_is_aligned(self, tmp_path):
        """minimap2 aligns the bundled aaftf/data/mito-start-cob.fasta to the assembly."""
        from importlib.resources import files

        from aaftf.mito import _orient_to_start

        fasta_in = tmp_path / "in.fa"
        fasta_in.write_text(">mt\nACGT\n")
        seen = []

        def _fake_paf_hits(cmd):
            seen.append(Path(cmd[-1]).read_text())
            return iter([])

        with patch("aaftf.mito.paf_hits", side_effect=_fake_paf_hits):
            _orient_to_start(str(fasta_in), str(tmp_path / "out.fa"))
        bundled = (files("aaftf") / "data" / "mito-start-cob.fasta").read_text()
        assert seen == [bundled] and bundled.startswith(">COB1")


class TestSubsample:
    """--subsample N keeps N read pairs with reformat.sh before NOVOPlasty runs."""

    def _run(self, tmp_path, subsample, reformat_rc=0):
        import subprocess

        workdir = tmp_path / "wd"
        cmds = []

        def _fake_run_cmd(cmd, debug=False, **kwargs):
            cmds.append(cmd)
            return subprocess.CompletedProcess(cmd, reformat_rc)

        with patch("aaftf.mito.require_tools"), patch("aaftf.mito.estimate_read_length", return_value=150), patch("aaftf.mito.run_cmd", side_effect=_fake_run_cmd):
            with patch("aaftf.mito.subprocess.Popen") as popen:
                popen.return_value.communicate.return_value = (b"", b"")
                with pytest.raises(RuntimeError):  # no NOVOPlasty output faked (or reformat.sh fails)
                    run(read1=str(tmp_path / "R1.fq"), read2=str(tmp_path / "R2.fq"), workdir=str(workdir), subsample=subsample)
        config = workdir / "novo-config.txt"
        return cmds, config.read_text() if config.exists() else None

    def test_subsampled_pairs_go_to_novoplasty(self, tmp_path):
        cmds, config = self._run(tmp_path, 2_000_000)
        assert cmds[0][0] == "reformat.sh"
        assert "samplereadstarget=2000000" in cmds[0] and any(a.startswith("sampleseed=") for a in cmds[0])
        assert f"in={tmp_path / 'R1.fq'}" in cmds[0] and f"in2={tmp_path / 'R2.fq'}" in cmds[0]
        assert "subsample_1.fastq.gz" in config and "subsample_2.fastq.gz" in config
        assert str(tmp_path / "R1.fq") not in config

    def test_default_keeps_1_5m_pairs(self):
        import inspect

        assert inspect.signature(run).parameters["subsample"].default == 1_500_000

    def test_zero_uses_all_reads(self, tmp_path):
        cmds, config = self._run(tmp_path, 0)
        assert cmds == [] and str(tmp_path / "R1.fq") in config

    def test_reformat_failure_raises(self, tmp_path):
        cmds, config = self._run(tmp_path, 1000, reformat_rc=1)
        assert len(cmds) == 1 and config is None  # stopped before writing the NOVOPlasty config

    def test_negative_subsample_raises(self, tmp_path):
        with pytest.raises(ValueError, match="--subsample"):
            run(read1="R1.fq", read2="R2.fq", workdir=str(tmp_path / "wd"), subsample=-1)


class TestNextCommandHint:
    """mito suggests filtering the same reads with the mitochondrial genome as -s/--screen_local."""

    def _run(self, tmp_path, caplog, **kwargs):
        import logging

        def _fake_popen(cmd, cwd, **kw):
            Path(cwd, "Circularized_assembly_1_x.fasta").write_text(">x\nACGT\n")
            return type("P", (), {"communicate": lambda self: None})()

        out = tmp_path / "strain.mito.fasta"
        with patch("aaftf.mito.require_tools"), patch("aaftf.mito.estimate_read_length", return_value=150):
            with patch("aaftf.mito.subprocess.Popen", side_effect=_fake_popen), patch("aaftf.mito._orient_to_start"):
                with patch("aaftf.mito._subsample_pairs", return_value=("wd/subsample_1.fastq.gz", "wd/subsample_2.fastq.gz")):
                    with caplog.at_level(logging.INFO, logger="aaftf"):
                        run(read1="trim/strain_1P.fastq.gz", read2="trim/strain_2P.fastq.gz", out=str(out), workdir=str(tmp_path / "wd"), **kwargs)
        return caplog.text, out

    def test_suggests_filter_with_mito_as_screen_local(self, tmp_path, caplog):
        text, out = self._run(tmp_path, caplog)
        assert f"AAFTF filter -1 trim/strain_1P.fastq.gz -2 trim/strain_2P.fastq.gz -s {out} -o strain" in text

    def test_hint_uses_original_reads_not_the_subsample(self, tmp_path, caplog):
        text, _ = self._run(tmp_path, caplog, subsample=1000)
        assert "-1 trim/strain_1P.fastq.gz" in text and "subsample_1" not in text.split("Your next command")[-1]

    def test_pipe_hides_hint(self, tmp_path, caplog):
        text, _ = self._run(tmp_path, caplog, subsample=0, pipe=True)
        assert "Your next command might be:" not in text
