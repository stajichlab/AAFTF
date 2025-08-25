"""fix tbl file after fcs_trimming.

This will attempt to adjust the coordinates for features in tbl file
after the contigs are modified.
"""
import csv
import re
import sys

from AAFTF.utility import status


def parse_tbl(tbl_file_handle):
    """Parse ncbi .tbl format file into a dictionary of lists."""
    features = {}
    sequence_name = None
    for line in tbl_file_handle:

        line = line.rstrip()
        if not line or line.startswith("#"):
            continue
        if line.startswith('>Feature'):
            m = re.match(r'>Feature\s+(\S+)', line)
            if m:
                sequence_name = m.group(1)
                features[sequence_name] = []
            else:
                status(f'Unexpected line {line} in tbl file')
                continue
        m = re.match(r'^([<>]?\d+)\t([<>]?\d+)(\t(\S+))?', line)
        if m:
            start = m.group(1)
            end = m.group(2)
            fname = m.group(4)
            if fname is not None:
                features[sequence_name].append([start, end, fname])
            else:
                features[sequence_name].append([start, end])
            continue
        m = re.match('^\t{3}(.+)', line)
        if m:
            features[sequence_name].append(['', '', '', m.group(1)])

    return features


def parse_adjustments(adj_file_handle):
    """Parse the adjustments listed from NCBI FCS trimming."""
    adjustments = {}

    ok = False
    for line in adj_file_handle:
        if line.startswith('#accession'):
            # we are in a proper format file
            ok = True
            continue
        if not ok:
            print(f'Skipping line: {line.strip()}', file=sys.stderr)
            continue
        break
    csvreader = csv.reader(adj_file_handle, delimiter="\t")
    for row in csvreader:
        if len(row) < 4:
            status(f'Skipping line: {row}')
            continue
        seqid = row[0]
        try:
            original_length = int(row[1])
            action = row[2]
            for range in row[3].split(','):
                m = re.match(r'(\d+)\.\.(\d+)', range)
                if m:
                    trim_start, trim_end = m.groups()
                    if seqid not in adjustments:
                        adjustments[seqid] = []
                    adjustments[seqid].append([int(original_length),
                                               action,
                                               int(trim_start), int(trim_end)])
        except ValueError:
            status(f'Skipping line: {row}')
            continue
    return adjustments


def fix_tbl(tbl_fh, adjustment_fh, output_handle):
    """Read table and apply adjustments to feature coordinates."""
    features = parse_tbl(tbl_fh)
    adjustments = parse_adjustments(adjustment_fh)
    for seqid, feats in features.items():
        print(f'>Feature {seqid}', file=output_handle)
        adj = {}
        if seqid in adjustments:
            for adjustment in adjustments[seqid]:
                original_length, action, trim_start, trim_end = adjustment
                if action == "ACTION_TRIM":
                    if trim_start == 1:
                        adj['trim_left'] = trim_end
                    elif trim_end == original_length:
                        adj['trim_right'] = trim_start
                    else:
                        status(f'Cannot trim effectively the adjustment is internal to the contig {seqid}:{trim_start}..{trim_end} in len={original_length}')
        for feature in feats:
            if feature[0] is not None and feature[0] != "":
                fstart = feature[0]
                fend = feature[1]
                if 'trim_left' in adj:
                    modstart = ""
                    modend = ""
                    m = re.match(r'([<>])', fstart)  # compile this for speed?
                    if m:
                        modstart = m.group(1)
                    m = re.match(r'([<>])', fend)  # compile this for speed?
                    if m:
                        modend = m.group(1)
                    fstart = int(fstart.lstrip('<>'))
                    fend = int(fend.lstrip('<>'))
                    fstart -= adj["trim_left"]
                    if fstart < 1:
                        fstart = 1
                    fend -= adj["trim_left"]
                    if fend < 1:
                        fend = 1
                    fstart = f'{modstart}{fstart}'
                    fend = f'{modend}{fend}'
                    feature[0] = fstart
                    feature[1] = fend
                if 'trim_right' in adj:
                    fstart = int(feature[0].lstrip('<>'))
                    fend = int(feature[1].lstrip('<>'))
                    if (fstart >= adj['trim_right'] or fend >= adj['trim_right']):
                        fend = adj['trim_right']
                        # status(f'Feature at {seqid}:{feature[0]}..{feature[1]} overlaps with right trim {adj["trim_right"]}')
            print("\t".join(feature), file=output_handle)


def run(parser, args):
    """Run the fix_tbl subcommand."""
    fix_tbl(args.table, args.report, args.output)
