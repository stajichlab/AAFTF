"""Run NCBI routines to identify contaminant and vectior contigs.

The contaminants are presumably sequences that were not screened out
in the filter step.

This uses NCBI fcs tool for screening which relies on a singularity
engine installed

The default libraries for screening are located in resources.py
and include common Euk, Prok, and MITO contaminants.
"""

import os
import shutil
import urllib.request
import uuid
from subprocess import call

from AAFTF.resources import FCSADAPTOR
from AAFTF.utility import SafeRemove, printCMD, status


def run(parser, args):
    """Perform vector trimming via the fcs screening tool."""
    DB = args.AAFTF_DB
    if args.AAFTF_DB:
        DB = args.AAFTF_DB
    elif "AAFTF_DB" in os.environ:
        DB = os.environ["AAFTF_DB"]
    else:
        DB = 'AAFTF_DB'
        if not os.path.exists(DB):
            os.mkdir(DB)
        status("No AAFTF_DB environ variable please see setup instructions. Creating locally.")

    containerengine = args.container_engine
    infilename = os.path.basename(os.path.realpath(args.infile))
    image = args.image
    tax = "--euk"
    if args.prok:
        tax = "--prok"
    custom_workdir = 1
    if not args.workdir:
        custom_workdir = 0
        args.workdir = 'aaftf-fcsscreen_'+str(uuid.uuid4())[:8]

    if not os.path.exists(args.workdir):
        os.mkdir(args.workdir)

    fcsexe = args.fcs_script
    if fcsexe is None:
        fcsexe = shutil.which('run_fcsadaptor.sh')
    if fcsexe is None:
        fcsexe = os.path.join(os.path.join(DB, 'run_fcsadaptor.sh'))
        #  This will help download the fcs-adaptor shell script rather than re-implementing it here
        if not os.path.exists(fcsexe):
            url = os.path.join(FCSADAPTOR['EXEURL'] % (FCSADAPTOR['VERSION']))
            if args.debug:
                status(f'url {url} download to {fcsexe}')
            urllib.request.urlretrieve(url, fcsexe)
            os.chmod(fcsexe, 0o755)
    if image is None:
        image = os.path.join(DB, FCSADAPTOR['SIFLOCAL'] % (FCSADAPTOR['VERSION']))
        if not os.path.exists(image):
            url = os.path.join(FCSADAPTOR['SIFURL'], FCSADAPTOR['VERSION'], FCSADAPTOR['SIF'])
            if args.debug:
                status(f'url {url} download to {image}')
            urllib.request.urlretrieve(url, image)

    cmd = [fcsexe, '--fasta-input', args.infile, '--output-dir', args.workdir, tax]

    if containerengine == "singularity":
        cmd += ['--container-engine', 'singularity', '--image', image]
    printCMD(cmd)
    DEVNULL = open(os.devnull, 'w')
    try:
        if args.debug:
            call(cmd)
        else:
            call(cmd, stderr=DEVNULL)
        DEVNULL.close()
    except NameError:
        print(f"error in calling executable {cmd}")

    cleanresult = os.path.join(args.workdir, 'cleaned_sequences', infilename)
    if args.debug:
        status(f'copy from: {cleanresult} -> {args.outfile}')
    os.rename(cleanresult, args.outfile)
    fcsreport = os.path.join(args.workdir, 'fcs_adaptor_report.txt')
    with open(fcsreport) as fh:
        status('FCS report:')
        for line in fh:
            print(line, end="")
    # make a copy of the report to show
    os.rename(fcsreport, args.outfile + ".fcs_adaptor_report.txt")
    # cleanup after running
    if not args.debug and not custom_workdir:
        SafeRemove(args.workdir)
