import os
import subprocess
import gzip
from Bio.SeqIO.FastaIO import SimpleFastaParser


def which_path(file_name):
    for path in os.environ["PATH"].split(os.pathsep):
        full_path = os.path.join(path, file_name)
        if os.path.exists(full_path) and os.access(full_path, os.X_OK):
            return full_path
    return None

def line_count(fname):
    with open(fname) as f:
        i = -1
        for i, l in enumerate(f):
            pass
    return i + 1

def countfasta(input):
    count = 0
    with open(input, 'rU') as f:
        for line in f:
            if line.startswith (">"):
                count += 1
    return count

def fastastats(input):
    count = 0
    length = 0
    with open(input, 'rU') as f:
    	for Header, Seq in SimpleFastaParser(f):
    		count += 1
    		length += len(Seq)
    return count,length

def countfastq(input):
    lines = sum(1 for line in open(input))
    count = int(lines) // 4
    return count
    
def softwrap(string, every=80):
    lines = []
    for i in range(0, len(string), every):
        lines.append(string[i:i+every])
    return '\n'.join(lines)
    
def bam_read_count(bamfile):
    cmd = ['samtools', 'idxstats', bamfile]
    mapped = 0
    unmapped = 0
    for line in execute(cmd, '.'):
        rname, rlen, nm, nu = line.rstrip().split()
        mapped += int(nm)
        unmapped += int(nu)
    return (mapped, unmapped)

def calcN50(lengths):
    lengths.sort()
    total_len = sum(lengths)
    n50 = 0
    cumulsum = 0
    i = 1
    for n in reversed(lengths):
        cumulsum += n
        if n50 == 0 and cumulsum >= total_len * 0.5:
            n50 = n
    return n50

def printCMD(cmd, num):
	return " \\\n\t\t".join([" ".join(cmd[i:i+num]) for i in range(0,len(cmd),num)])

#from https://stackoverflow.com/questions/4417546/constantly-print-subprocess-output-while-process-is-running
def execute(cmd, dir):
    DEVNULL = open(os.devnull, 'w')
    popen = subprocess.Popen(cmd, cwd=dir, stdout=subprocess.PIPE, universal_newlines=True, stderr=DEVNULL)
    for stdout_line in iter(popen.stdout.readline, ""):
        yield stdout_line 
    popen.stdout.close()
    return_code = popen.wait()
    if return_code:
        raise subprocess.CalledProcessError(return_code, cmd)
