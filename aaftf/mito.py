"""Run the Mitochondria assembly tool NOVOPlasty."""

import logging
import os
import subprocess
import sys
import uuid
from importlib.resources import as_file, files
from pathlib import Path
from typing import Any

from Bio.SeqIO.FastaIO import SimpleFastaParser

from aaftf.utility import COMPLEMENT, basename_from_reads, cleanup_workdir, close_step_log, estimate_read_length, make_workdir, paf_hits, print_cmd, require_tools, run_cmd, write_fasta

__all__ = ["run"]


logger = logging.getLogger(__name__)

_PACKAGE_DATA = files("aaftf") / "data"
# fixed reformat.sh seed so --subsample picks the same read pairs on every run
_SUBSAMPLE_SEED = 42


def run(
    read1: str,
    read2: str,
    out: str = "mito.fasta",
    workdir: str | None = None,
    minlen: int = 10000,
    maxlen: int = 100000,
    seed: str | None = None,
    subsample: int = 1_500_000,
    memory: int = 8,
    debug: bool = False,
    pipe: bool = False,
    **kwargs: Any,
) -> None:
    """Assemble a mitochondrial genome from paired reads with NOVOPlasty.

    Two sequences steer the run:

    * the **seed** is where NOVOPlasty starts assembling; it extends the seed with matching
      reads until the genome is complete (bundled default: an *A. nidulans* cob fragment,
      ``aaftf/data/mito-seed.fasta``);
    * the **start gene** only matters afterwards: a circular assembly is rotated to begin at the
      cob gene, found by aligning a bundled consensus of fungal cob genes
      (``aaftf/data/mito-start-cob.fasta``).

    Mitochondrial reads are usually at far higher coverage than nuclear ones, so by default the
    input is capped at ``subsample`` randomly chosen read pairs (the same pairs on every run) to
    save time and memory. NOVOPlasty may also subsample on its own to stay within ``memory``.

    A circularized assembly is written as a single rotated ``mt`` record; otherwise the contigs
    are written as ``contig_N``. The suggested next command is ``filter`` on the same reads with the
    mitochondrial genome added via ``-s/--screen_local``, so mitochondrial reads are screened out
    before the nuclear assembly.

    Args:
        read1: Forward reads FASTQ.
        read2: Reverse reads FASTQ.
        out: Output FASTA path.
        workdir: Working directory; a temporary one is created when None.
        minlen: Minimum expected genome size (NOVOPlasty genome range).
        maxlen: Maximum expected genome size (NOVOPlasty genome range).
        seed: FASTA of the NOVOPlasty seed; None uses the bundled *A. nidulans* cob fragment.
        subsample: Keep only this many randomly chosen read pairs (BBTools ``reformat.sh``);
            0 uses all reads.
        memory: Max memory in GB for NOVOPlasty (and ``reformat.sh``).
        debug: Keep the work directory.
        pipe: Suppress the "next command" hint; set by the ``pipeline_*`` subcommands (not a CLI option).
        **kwargs: Other parsed CLI attributes (``command``, ``func``, ``quiet``, ...); ignored.

    Raises:
        ValueError: If ``subsample`` is negative.
        RuntimeError: If subsampling fails or NOVOPlasty produces no assembly.
    """
    if subsample < 0:
        raise ValueError(f"--subsample must be a number of read pairs, or 0 to use all reads; got {subsample}")
    require_tools(["NOVOPlasty.pl", "minimap2"] + (["reformat.sh"] if subsample else []))
    # first we need to generate working directory
    unique_id = str(uuid.uuid4())[:8]
    workdir, custom_workdir = make_workdir(workdir, "mito")
    try:
        input_reads = (read1, read2)  # the hint suggests filtering these, not the subsample
        if subsample:
            read1, read2 = _subsample_pairs(read1, read2, subsample, workdir, memory, debug)

        # now estimate read lengths of FASTQ
        read_len = estimate_read_length(read1)

        # NOVOPlasty seed: --seed, else the bundled default (copied into the work directory so
        # NOVOPlasty gets a real file path even from a zipped install)
        if seed:
            seed_fasta = str(Path(seed).resolve())
        else:
            seed_fasta = str(Path(workdir, "mito-seed.fasta").resolve())
            Path(seed_fasta).write_bytes((_PACKAGE_DATA / "mito-seed.fasta").read_bytes())

        # write the NOVOPlasty config from the bundled template
        placeholders = {
            "<PROJECT>": unique_id,
            "<MINLEN>": str(minlen),
            "<MAXLEN>": str(maxlen),
            "<MAXMEM>": str(memory),
            "<SEED>": seed_fasta,
            "<READLEN>": str(read_len),
            "<FORWARD>": str(Path(read1).resolve()),
            "<REVERSE>": str(Path(read2).resolve()),
        }
        config_text = (_PACKAGE_DATA / "novoplasty-config.txt").read_text()
        for placeholder, value in placeholders.items():
            config_text = config_text.replace(placeholder, value)
        Path(workdir, "novo-config.txt").write_text(config_text)

        # now we can finally run NOVOplasty.pl
        logger.info("De novo assembling mitochondrial genome using NOVOplasty")
        cmd = ["NOVOPlasty.pl", "-c", "novo-config.txt"]
        print_cmd(cmd)
        novolog = str(Path(workdir, "novoplasty.log"))
        with open(novolog, "w") as logfile:
            p1 = subprocess.Popen(cmd, cwd=workdir, stdout=logfile, stderr=logfile)
            p1.communicate()

        # now parse the results, preferring a circular assembly over partial ones
        outputs = sorted(os.listdir(workdir))
        draft_mito = None
        for output_prefix in ("Circularized_assembly_", "Contigs_1_", "Uncircularized_assemblies_"):
            draft_mito = next((str(Path(workdir, f)) for f in outputs if f.startswith(output_prefix)), None)
            if draft_mito:
                break
        if draft_mito is None:
            raise RuntimeError(f"NOVOplasty did not produce an assembly - check {novolog}")
        circular = Path(draft_mito).name.startswith("Circularized_assembly_")
        if circular:
            logger.info("NOVOplasty assembled complete circular genome")
            logger.info("Rotating assembly to start at the cytochrome b (cob) gene")
            _orient_to_start(draft_mito, out)
        else:
            num_contigs = 0
            contig_length = 0
            with open(out, "w") as outfile:
                with open(draft_mito) as infile:
                    for title, seq in SimpleFastaParser(infile):
                        num_contigs += 1
                        contig_length += len(seq)
                        write_fasta(outfile, f"contig_{num_contigs}", seq)
            # weird formatting here for PEP8
            logger.info(f"NOVOplasty assembled {num_contigs} contigs consisting of {contig_length:,} bp," + "but was unable to circularize genome")

        logger.info(f"AAFTF mito complete: {out}")
        cleanup_workdir(workdir, debug, custom_workdir)
        if not pipe:
            read1, read2 = input_reads
            logger.info(f"Your next command might be:\nAAFTF filter -1 {read1} -2 {read2} -s {out} -o {basename_from_reads(read1)}")
    finally:
        close_step_log()


def _subsample_pairs(read1: str, read2: str, pairs: int, workdir: str, memory: int, debug: bool) -> tuple[str, str]:
    """Randomly keep ``pairs`` read pairs from ``read1``/``read2`` with BBTools ``reformat.sh``.

    Both files are sampled together so mates stay paired, with a fixed seed so reruns keep the
    same pairs. If there are fewer than ``pairs`` pairs, all of them are kept.

    Args:
        read1: Forward reads FASTQ.
        read2: Reverse reads FASTQ.
        pairs: Number of read pairs to keep.
        workdir: Directory the subsampled FASTQ files are written to.
        memory: Java heap size for ``reformat.sh``, in GB.
        debug: Show ``reformat.sh`` output.

    Returns:
        The subsampled forward and reverse FASTQ paths.

    Raises:
        RuntimeError: If ``reformat.sh`` fails.
    """
    sub_left = str(Path(workdir, "subsample_1.fastq.gz"))
    sub_right = str(Path(workdir, "subsample_2.fastq.gz"))
    cmd = [
        "reformat.sh",
        f"-Xmx{memory}g",
        f"in={read1}",
        f"in2={read2}",
        f"out={sub_left}",
        f"out2={sub_right}",
        f"samplereadstarget={pairs}",
        f"sampleseed={_SUBSAMPLE_SEED}",
        "overwrite=t",
    ]
    logger.info(f"Subsampling to {pairs:,} read pairs")
    if run_cmd(cmd, debug).returncode != 0:
        raise RuntimeError(f"reformat.sh failed to subsample {read1} / {read2}")
    return sub_left, sub_right


def _orient_to_start(fasta_in: str, fasta_out: str) -> None:
    """Rotate a circular MT assembly to begin at the cob gene and write it as ``mt``.

    The bundled consensus of fungal cob genes (``aaftf/data/mito-start-cob.fasta``) is aligned
    to the assembly with minimap2; with exactly one in-range hit
    the sequence is rotated (and reverse complemented for a minus-strand hit). Otherwise the
    last sequence of ``fasta_in`` is written unrotated and an error is logged.

    Args:
        fasta_in: FASTA of the circular assembly.
        fasta_out: Output FASTA path.
    """
    # load sequence into dictionary
    initial_seq = ""
    # header = ''
    with open(fasta_in) as infile:
        for title, seq in SimpleFastaParser(infile):
            initial_seq = seq
            # header = title

    # align the bundled spoa consensus of fungal cob genes; as_file gives it a real path (a
    # temporary copy only when the package is installed zipped)
    with as_file(_PACKAGE_DATA / "mito-start-cob.fasta") as cob_fasta:
        alignments = list(paf_hits(["minimap2", "-x", "map-ont", "-c", fasta_in, str(cob_fasta)]))
    if len(alignments) == 1:
        hit = alignments[0]
        ref_strand = hit.strand
        ref_offset = hit.query_start
        if ref_strand == "-":
            ref_start = hit.target_end + ref_offset
        else:
            ref_start = hit.target_start - ref_offset
        if ref_start == len(initial_seq):
            ref_start = 0
        if ref_start < 0 or ref_start > len(initial_seq):
            # A partial alignment of the seed to the contig edge can push
            # this offset out of range; Python's negative-index slicing
            # would silently wrap and produce a bogus rotation instead of
            # erroring, so treat it the same as a failed rotation.
            logger.error(f"unable to rotate because computed rotation offset {ref_start} is out of range for sequence of length {len(initial_seq)}")
            with open(fasta_out, "w") as outfile:
                write_fasta(outfile, "mt", initial_seq)
            return
        rotated = initial_seq[ref_start:] + initial_seq[:ref_start]
        if ref_strand == "-":
            rotated = rotated.translate(COMPLEMENT)[::-1]  # reverse complement
        with open(fasta_out, "w") as outfile:
            write_fasta(outfile, "mt", rotated)
    elif len(alignments) == 0:
        logger.error("unable to rotate because the cob gene was not found in the assembly")
        with open(fasta_out, "w") as outfile:
            write_fasta(outfile, "mt", initial_seq)
    elif len(alignments) > 1:
        logger.error("unable to rotate because found multiple alignments")
        for x in alignments:
            sys.stderr.write(f"{x}\n")
        with open(fasta_out, "w") as outfile:
            write_fasta(outfile, "mt", initial_seq)
