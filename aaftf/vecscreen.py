"""Identify and remove contaminant, vector and mitochondrial contigs from an assembly.

The contaminants are presumably sequences that were not screened out
in the filter step. Contigs are screened with BLASTN against the UniVec
vector library and against eukaryotic, prokaryotic and mitochondrial
contaminant databases.

The default libraries for screening are located in resources.py
and include common Euk, Prok, and MITO contaminants.
"""

import csv
import logging
import os
from collections.abc import Iterator
from pathlib import Path
from subprocess import DEVNULL, call
from typing import Any

# biopython needed
from Bio import SeqIO

from aaftf.utility import cleanup_workdir, close_step_log, concat_files, make_workdir, next_step_name, print_cmd, require_databases, write_fasta

__all__ = ["BLAST_PERCENT_ID_CONTAM_MATCH", "BLAST_PERCENT_ID_MITO_MATCH", "run"]


BLAST_PERCENT_ID_CONTAM_MATCH = "90.0"

BLAST_PERCENT_ID_MITO_MATCH = "98.6"

# VecScreen matches
"""
Strong Match to Vector
(Expect 1 random match in 1,000,000 queries of length 350 kb.)
Terminal match with Score ≥ 24.
Internal match with Score ≥ 30.
Moderate Match to Vector
(Expect 1 random match in 1,000 queries of length 350 kb.)
Terminal match with Score 19 to 23.
Internal match with Score 25 to 29.
Weak Match to Vector
(Expect 1 random match in 40 queries of length 350 kb.)
Terminal match with Score 16 to 18.
Internal match with Score 23 to 24.
Segment of Suspect Origin
Any segment of fewer than 50 bases between two vector matches
    or between a match and an end.
"""

logger = logging.getLogger(__name__)

# BLAST database name used by the screens -> `AAFTF database` database it is built from
_CONTAM_BLAST_DBS = {"UniVec": "univec", "CONTAM_EUKS": "euks", "CONTAM_PROKS": "proks", "MITO": "mitodb"}


def run(
    infile: str,
    outfile: str,
    workdir: str | None = None,
    cpus: int = 1,
    percent_id: str | None = None,
    stringency: str = "high",
    debug: bool = False,
    pipe: bool = False,
    **kwargs: Any,
) -> None:
    """Screen an assembly for contaminants and vectors with BLASTN.

    Pipeline: build the contamination BLAST databases, screen out Euk/Prok
    contaminants, screen for mitochondrial contigs, run repeated VecScreen
    (UniVec) rounds to trim/split vector hits, then write the cleaned
    assembly and a separate mitochondrial-contigs FASTA.

    Args:
        infile: Input assembly FASTA.
        outfile: Output cleaned FASTA; only its basename is used, so it is written to the
            current directory.
        workdir: Working directory; a temporary one is created if None.
        cpus: Number of BLAST threads.
        percent_id: BLASTN ``-perc_identity`` cutoff (as a string) for the Euk/Prok screen;
            defaults to ``BLAST_PERCENT_ID_CONTAM_MATCH``.
        stringency: ``"high"`` keeps moderate and strong vector hits; anything else keeps
            only strong hits.
        debug: Keep the working directory when True.
        pipe: Suppress the "next command" hint; set by the ``pipeline_*`` subcommands (not a CLI option).
        **kwargs: Other parsed CLI attributes (``command``, ``func``, ``quiet``); ignored.
    """
    workdir, custom_workdir = make_workdir(workdir, "vecscreen")
    try:
        percentid_cutoff = percent_id or BLAST_PERCENT_ID_CONTAM_MATCH

        # final_outfile/outdir/prefix are derived from the user's --outfile once,
        # up front, so nothing later in the pipeline can accidentally clobber them
        # (each stage below only sees its own function-local variables).
        final_outfile = outfile
        outdir = str(Path(outfile).parent)
        outname = Path(outfile).name
        if ".f" in outname:
            prefix = outname.rsplit(".f", 1)[0]
        else:
            prefix = str(os.getpid())
        _build_contam_databases(workdir)

        contigs_to_remove: dict[str, tuple[str, str, float]] = {}
        regions_to_trim = _screen_euk_prok_contamination(infile, workdir, prefix, cpus, percentid_cutoff)
        euk_cleaned = _write_euk_cleaned(infile, regions_to_trim, workdir, prefix)
        mito_hits = _screen_mitochondria(euk_cleaned, workdir, prefix, cpus, contigs_to_remove)
        outfile_vec = _run_vecscreen_rounds(euk_cleaned, workdir, prefix, cpus, stringency, contigs_to_remove)

        logger.info(f"{len(contigs_to_remove):,} contigs will be removed:")
        for k, v in sorted(contigs_to_remove.items()):
            print(f"\t{k} --> dbhit={v[0]}; hit={v[1]}; pident={v[2]}")

        # this could instead use the outfile and strip
        # .fasta/fsa/fna and add mito on it I suppose, but assumes
        # a bit about the naming structure
        mitochondria = str(Path(outdir, prefix + ".mitochondria.fasta"))
        _write_final_outputs(outfile_vec, contigs_to_remove, mito_hits, final_outfile, mitochondria)

        next_out = next_step_name(final_outfile, ".taxpurge.fasta")
        if not pipe:
            logger.info("Your next command might be:\n" + "AAFTF taxpurge -i {:} -o {:} --phylum {:}".format(final_outfile, next_out, "Ascomycota"))

        cleanup_workdir(workdir, debug, custom_workdir)
    finally:
        close_step_log()


def _build_contam_databases(workdir: str) -> None:
    """Build a BLAST nucleotide DB in ``workdir`` for each contamination screen.

    Args:
        workdir: Directory in which the combined FASTA files and BLAST DBs are written.
    """
    logger.info("Building BLAST databases for contamination screen.")
    sources = require_databases(list(_CONTAM_BLAST_DBS.values()))
    for blast_name, source in zip(_CONTAM_BLAST_DBS, sources):
        combined_fasta = str(Path(workdir, f"{blast_name}.fasta"))
        concat_files([source], combined_fasta)
        _make_blastdb("nucl", combined_fasta, str(Path(workdir, blast_name)))


def _make_blastdb(type: str, file: str, name: str) -> None:
    """Run ``makeblastdb`` unless an up-to-date index already exists.

    Args:
        type: BLAST DB type, ``"nucl"`` or ``"prot"``.
        file: Input FASTA file.
        name: Output database name (path prefix).
    """
    idxfile = name
    if type == "nucl":
        idxfile += ".nin"
    else:
        idxfile += ".pin"
    idxexists = Path(idxfile).exists()
    if not idxexists or Path(idxfile).stat().st_ctime < Path(file).stat().st_ctime:
        cmd = ["makeblastdb", "-dbtype", type, "-in", file, "-out", name]
        print_cmd(cmd)
        call(cmd, stdout=DEVNULL, stderr=DEVNULL)


def _screen_euk_prok_contamination(infile: str, workdir: str, prefix: str, cpus: int, percentid_cutoff: str) -> dict[str, list[tuple[int, int, str, str, float]]]:
    """BLASTN ``infile`` against the Euk/Prok contamination DBs and collect qualifying hits.

    A hit qualifies at >=98% identity over >=50 bp, >=94% over >=100 bp, or >=90% over >=200 bp.

    Args:
        infile: Query assembly FASTA.
        workdir: Directory holding the BLAST DBs and receiving the reports.
        prefix: Prefix for the BLAST report file names.
        cpus: Number of BLAST threads.
        percentid_cutoff: BLASTN ``-perc_identity`` value.

    Returns:
        Mapping of contig id to a list of ``(start, end, db_name, hit_id, pident)`` tuples.
    """
    # qaccver saccver pident length mismatch gapopen qstart qend
    # sstart send evalue bitscore
    regions_to_trim: dict[str, list[tuple[int, int, str, str, float]]] = {}
    for contam in ["CONTAM_EUKS", "CONTAM_PROKS"]:
        logger.info(f"{contam} Contamination Screen")
        for row in _run_blastn_screen(infile, workdir, prefix, contam, cpus, percentid_cutoff):
            if (float(row[2]) >= 98.0 and int(row[3]) >= 50) or (float(row[2]) >= 94.0 and int(row[3]) >= 100) or (float(row[2]) >= 90.0 and int(row[3]) >= 200):
                start, end = sorted([int(row[6]), int(row[7])])
                regions_to_trim.setdefault(row[0], []).append((start, end, contam, row[1], float(row[2])))
        logger.info(f"{contam} screening finished")
    return regions_to_trim


def _run_blastn_screen(query: str, workdir: str, prefix: str, dbname: str, cpus: int, percent_identity: str) -> list[list[str]]:
    """Run blastn (tab-6 output) for ``query`` against a DB built in ``workdir``.

    Args:
        query: Query FASTA file.
        workdir: Directory holding the DB and receiving the report.
        prefix: Prefix used in the report file name.
        dbname: Name of the BLAST DB inside ``workdir``.
        cpus: Number of BLAST threads.
        percent_identity: BLASTN ``-perc_identity`` value.

    Returns:
        The report rows, each a list of column strings.
    """
    blastreport = str(Path(workdir, f"{dbname}.{prefix}.blastn"))
    blastnargs = ["blastn", "-query", query, "-db", str(Path(workdir, dbname)), "-num_threads", str(cpus), "-dust", "yes", "-soft_masking", "true", "-perc_identity", percent_identity, "-lcase_masking", "-outfmt", "6", "-out", blastreport]
    print_cmd(blastnargs)
    call(blastnargs)
    with open(blastreport) as report:
        return list(csv.reader(report, delimiter="\t"))


def _write_euk_cleaned(infile: str, regions_to_trim: dict[str, list[tuple[int, int, str, str, float]]], workdir: str, prefix: str) -> str:
    """Split out Euk/Prok-contaminated regions found by ``_screen_euk_prok_contamination``.

    Args:
        infile: Input assembly FASTA.
        regions_to_trim: Contaminated regions per contig id.
        workdir: Directory for the cleaned FASTA.
        prefix: Prefix for the cleaned FASTA file name.

    Returns:
        Path to the cleaned FASTA, or ``infile`` unchanged if there was nothing to trim.
    """
    if not regions_to_trim:
        return infile

    euk_cleaned = str(Path(workdir, f"{prefix}.euk-prot_cleaned.fasta"))
    with open(euk_cleaned, "w") as cleanout, open(infile) as fastain:
        for record in SeqIO.parse(fastain, "fasta"):
            if record.id not in regions_to_trim:
                write_fasta(cleanout, record.id, str(record.seq))
            else:
                contig_seq = str(record.seq)
                regions = regions_to_trim[record.id]
                logger.info(f"Splitting {record.id} for contamination: {regions}")
                lastpos = 0
                for i, x in enumerate(regions):
                    # x[0] is a 1-based BLAST start; subtract one so the
                    # slice end doesn't retain the first contaminant base.
                    new_seq = contig_seq[lastpos : x[0] - 1]
                    lastpos = x[1]
                    write_fasta(cleanout, f"split{i}_{record.id}", new_seq)
                    if i == len(regions) - 1:
                        new_seq = contig_seq[x[1] :]
                        write_fasta(cleanout, f"split{i + 1}_{record.id}", new_seq)
    return euk_cleaned


def _screen_mitochondria(euk_cleaned: str, workdir: str, prefix: str, cpus: int, contigs_to_remove: dict[str, tuple[str, str, float]]) -> list[str]:
    """BLASTN against the MITO DB and flag hits of >=120 bp for removal.

    Args:
        euk_cleaned: Query FASTA (output of the Euk/Prok screen).
        workdir: Directory holding the DB and receiving the report.
        prefix: Prefix for the report file name.
        cpus: Number of BLAST threads.
        contigs_to_remove: Mapping of contig id to ``(screen, hit_id, pident)``; mutated in place.

    Returns:
        Contig ids identified as mitochondrial.
    """
    logger.info("Mitochondria Contamination Screen")
    mito_hits = []
    for row in _run_blastn_screen(euk_cleaned, workdir, prefix, "MITO", cpus, BLAST_PERCENT_ID_MITO_MATCH):
        if int(row[3]) >= 120:
            contigs_to_remove[row[0]] = ("MitoScreen", row[1], float(row[2]))
            mito_hits.append(row[0])
    logger.info("Mito screening finished.")
    return mito_hits


def _run_vecscreen_rounds(euk_cleaned: str, workdir: str, prefix: str, cpus: int, stringency: str, contigs_to_remove: dict[str, tuple[str, str, float]]) -> str:
    """Repeatedly BLASTN against UniVec, trimming/splitting vector hits each round until none remain.

    An existing round report in ``workdir`` is reused rather than recomputed.

    Args:
        euk_cleaned: Starting query FASTA.
        workdir: Directory holding the UniVec DB and round outputs.
        prefix: Prefix for the per-round file names.
        cpus: Number of BLAST threads.
        stringency: Vector-hit stringency passed to ``_parse_clean_blastn``.
        contigs_to_remove: Contigs already flagged for removal; their hits are skipped.

    Returns:
        Path to the final, fully vector-cleaned FASTA.
    """
    logger.info("Starting VecScreen, will remove terminal matches and split internal matches")
    rnd = 0
    count = 1
    cleanfile = euk_cleaned
    while count > 0:
        filepref = f"{prefix}.r{rnd}"
        report = str(Path(workdir, f"{filepref}.vecscreen.tab"))
        if not Path(report).exists():
            cmd = [
                "blastn",
                "-task",
                "blastn",
                "-reward",
                "1",
                "-penalty",
                "-5",
                "-gapopen",
                "3",
                "-gapextend",
                "3",
                "-dust",
                "yes",
                "-soft_masking",
                "true",
                "-evalue",
                "700",
                "-searchsp",
                "1750000000000",
                "-db",
                str(Path(workdir, "UniVec")),
                "-outfmt",
                "6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore score qlen",
                "-num_threads",
                str(cpus),
                "-query",
                euk_cleaned,
                "-out",
                report,
            ]
            print_cmd(cmd)
            call(cmd)
        logger.info(f"Parsing VecScreen round {rnd + 1}: {filepref} for {report}")

        (count, cleanfile) = _parse_clean_blastn(euk_cleaned, str(Path(workdir, filepref)), report, stringency, contigs_to_remove)
        logger.info(f"count is {count} cleanfile is {cleanfile}")
        if count > 0:  # vector matches remain; re-screen the newly trimmed/split sequences
            rnd += 1
            euk_cleaned = cleanfile

    return cleanfile


def _parse_clean_blastn(fastafile: str, prefix: str, blastn: str, stringent: str, contigs_to_remove: dict[str, tuple[str, str, float]] | None = None) -> tuple[int, str]:
    """Parse a VecScreen BLASTN report and write a vector-trimmed/split FASTA.

    Args:
        fastafile: FASTA that was screened.
        prefix: Path prefix; the output is written to ``<prefix>.clean.fsa``.
        blastn: VecScreen BLASTN tab report.
        stringent: ``"high"`` keeps moderate and strong hits, otherwise only strong ones.
        contigs_to_remove: Contigs whose hits are skipped; None means none.

    Returns:
        Tuple of (number of qualifying vector hits, cleaned FASTA path).
    """
    cleaned = prefix + ".clean.fsa"
    if contigs_to_remove is None:
        contigs_to_remove = {}
    vec_hits, found_vector_seq = _classify_vector_hits(blastn, stringent, contigs_to_remove)
    _write_trimmed_and_split(fastafile, vec_hits, cleaned)
    return found_vector_seq, cleaned


def _classify_vector_hits(blastn: str, stringent: str, contigs_to_remove: dict[str, tuple[str, str, float]]) -> tuple[dict[str, list[tuple[str, int, list[int], int, bool, str | None]]], int]:
    """Read a VecScreen-style BLASTN tab report and classify each surviving hit.

    Input rows have columns: qaccver saccver pident length mismatch gapopen
    qstart qend sstart send evalue bitscore score qlen

    Each hit is classified as terminal (within 25bp of either end of the
    query) or internal, then scored weak/moderate/strong per the VecScreen
    table (see module docstring):
    https://www.ncbi.nlm.nih.gov/tools/vecscreen/about/#Moderate
    (using ``score``, not ``bitscore``). Hits on a contig already in
    ``contigs_to_remove`` are skipped. Weak hits are always dropped;
    ``stringent == "high"`` additionally keeps moderate hits, otherwise only
    strong hits are kept.

    Args:
        blastn: VecScreen BLASTN tab report.
        stringent: ``"high"`` keeps moderate and strong hits, otherwise only strong ones.
        contigs_to_remove: Contigs whose hits are skipped.

    Returns:
        Tuple ``(vec_hits, found_vector_seq)`` where ``vec_hits`` maps contig id to a list of
        ``(hit_id, qlen, loc, score, terminal, position)`` and ``found_vector_seq`` is the
        number of qualifying hits.
    """
    vec_hits: dict[str, list[tuple[str, int, list[int], int, bool, str | None]]] = {}
    found_vector_seq = 0
    with open(blastn) as vectab:
        rdr = csv.reader(vectab, delimiter="\t")
        for row in rdr:
            (qaccver, saccver, pid, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore, score_str, qlen) = row
            if qaccver in contigs_to_remove:
                continue

            loc = sorted([int(qstart), int(qend)])
            terminal = False
            position = None
            if loc[0] <= 25:
                terminal = True
                position = "5"
            if (int(qlen) - loc[1]) <= 25:
                terminal = True
                position = "3"

            score = int(score_str)
            match_strength = 0  # weak=0, moderate=1, strong=2
            if terminal:
                if score >= 19:
                    match_strength = 1
                if score >= 24:
                    match_strength = 2
            else:
                if score >= 25:
                    match_strength = 1
                if score >= 30:
                    match_strength = 2

            if match_strength == 0:
                continue
            if stringent != "high" and match_strength < 2:
                continue

            found_vector_seq += 1
            vec_hits.setdefault(qaccver, []).append((saccver, int(qlen), loc, score, terminal, position))
    return vec_hits, found_vector_seq


# flake8: noqa: C901
def _write_trimmed_and_split(fastafile: str, vec_hits: dict[str, list[tuple[str, int, list[int], int, bool, str | None]]], cleaned: str) -> None:
    """Write ``cleaned`` FASTA, trimming terminal vector hits and splitting out internal ones.

    For each record with hits in ``vec_hits``: terminal 5' hits push the
    kept region's start forward, terminal 3' hits pull its end back, and
    internal hits carve the contig into multiple ``splitN_<id>`` pieces
    around the vector regions (any hits within 50bp of each other collapse
    into one cut, handled naturally across ``_run_vecscreen_rounds``'
    repeated rounds rather than here). Records/pieces shorter than 200bp are
    dropped.

    Args:
        fastafile: Input FASTA.
        vec_hits: Vector hits per contig id, from ``_classify_vector_hits``.
        cleaned: Output FASTA path.
    """
    with open(cleaned, "w") as output_handle, open(fastafile) as fastain:
        for record in SeqIO.parse(fastain, "fasta"):
            seq_str = str(record.seq)
            if record.id not in vec_hits:
                if len(record.seq) >= 200:
                    write_fasta(output_handle, record.id, seq_str)
                continue

            five_end = 0
            three_end = len(record.seq)
            internal_regions = []
            for hit_id, qlen, loc, score, terminal, position in vec_hits[record.id]:
                if terminal and position == "5":
                    five_end = max(five_end, loc[1])
                elif terminal and position == "3":
                    # loc[0] (qstart) is 1-based; the slice end must be
                    # one less so the first vector/contaminant base isn't kept.
                    three_end = min(three_end, loc[0] - 1)
                elif loc not in internal_regions:
                    internal_regions.append(loc)

            sorted_internals = sorted(internal_regions, key=lambda region: region[0])
            slicer = [five_end]
            for start, end in sorted_internals:
                # start (qstart) is 1-based and used as a slice end below, so
                # shift by one to avoid retaining the first hit base. end
                # (qend) is used as the following slice start, already correct.
                slicer += [start - 1, end]
            slicer.append(three_end)
            keep_regions = list(_group(slicer, 2))

            if len(keep_regions) < 2:
                logger.info(f"Terminal trimming {record.id} to {keep_regions}")
                new_seq = seq_str[keep_regions[0][0] : keep_regions[0][1]]
                if len(new_seq) >= 200:
                    write_fasta(output_handle, record.id, new_seq)
            else:
                logger.info(f"Splitting contig {record.id} into {keep_regions}")
                for num, (start, end) in enumerate(keep_regions):
                    new_seq = seq_str[start:end]
                    if len(new_seq) >= 200:
                        write_fasta(output_handle, f"split{num + 1}_{record.id}", new_seq)


def _group(lst: list[int], n: int) -> Iterator[tuple[int, ...]]:
    """Yield consecutive non-overlapping ``n``-sized tuples from ``lst``, dropping a short tail.

    Args:
        lst: Values to group.
        n: Group size.

    Yields:
        Tuples of ``n`` consecutive values.
    """
    for i in range(0, len(lst), n):
        val = lst[i : i + n]
        if len(val) == n:
            yield tuple(val)


def _write_final_outputs(outfile_vec: str, contigs_to_remove: dict[str, tuple[str, str, float]], mito_hits: list[str], outfile: str, mitochondria: str) -> None:
    """Split the vecscreen output into the cleaned assembly and a mitochondrial-contigs FASTA.

    Args:
        outfile_vec: Vector-cleaned FASTA.
        contigs_to_remove: Contigs excluded from the cleaned assembly.
        mito_hits: Removed contigs that are written to the mitochondrial FASTA.
        outfile: Output cleaned assembly FASTA.
        mitochondria: Output mitochondrial-contigs FASTA.
    """
    n_clean = n_mito = 0
    with open(outfile, "w") as oh, open(mitochondria, "w") as mh, open(outfile_vec) as vec_in:
        for record in SeqIO.parse(vec_in, "fasta"):
            if record.id not in contigs_to_remove:
                write_fasta(oh, record.description, str(record.seq))
                n_clean += 1
            elif record.id in mito_hits:
                write_fasta(mh, record.description, str(record.seq))
                n_mito += 1
    logger.info(f"Writing {n_clean:,} cleaned contigs to: {outfile}")
    logger.info(f"Writing {n_mito:,} mitochondrial contigs to: {mitochondria}")
