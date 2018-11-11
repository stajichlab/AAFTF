import os
import subprocess
from Bio.SeqIO.FastaIO import SimpleFastaParser
import shutil
import textwrap
import datetime

def getRAM():
    import resource
    import sys
    rusage_denom = 1024.
    if sys.platform == 'darwin':
        # ... it seems that in OSX the output is different units ...
        rusage_denom = rusage_denom * rusage_denom
    mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / rusage_denom
    return mem

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
    lines = sum(1 for line in zopen(input))
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

def calcN50(lengths, num=0.5):
    lengths.sort()
    total_len = sum(lengths)
    n50 = 0
    cumulsum = 0
    i = 1
    for n in reversed(lengths):
        cumulsum += n
        if n50 == 0 and cumulsum >= total_len * num:
            n50 = n
    return n50


def printCMD(cmd):
	stringcmd = '{:}'.format(' '.join(cmd))
	prefix = '\033[96mCMD:\033[00m '
	wrapper = textwrap.TextWrapper(initial_indent=prefix, width=80, subsequent_indent=' '*8, break_long_words=False)
	print(wrapper.fill(stringcmd))

def status(string):
	print('\033[92m[{:}]\033[00m {:}'.format(datetime.datetime.now().strftime('%b %d %I:%M %p'), string))

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

def Fzip_inplace(input, cpus):
    '''
    function to zip as fast as it can, pigz -> gzip
    '''
    if which_path('pigz'):
        cmd = ['pigz', '-f', '-p', str(cpus), input]
    else:
        cmd = ['gzip', '-f', input]
    try:
        runSubprocess(cmd, '.', log)
    except NameError:
        subprocess.call(cmd)

def SafeRemove(input):
    if os.path.isdir(input):
        shutil.rmtree(input)
    elif os.path.isfile(input):
        os.remove(input)
    else:
        return
        
#streaming parallel pigz open via https://github.com/DarkoVeberic/utl/blob/master/futile/futile.py 
def which(program):
    import os
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)
    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None

def open_pipe(command, mode='r', buff=1024*1024):
    import subprocess
    import signal
    if 'r' in mode:
        return subprocess.Popen(command, shell=True, bufsize=buff,
                                stdout=subprocess.PIPE,
                                preexec_fn=lambda: signal.signal(signal.SIGPIPE, signal.SIG_DFL)
                               ).stdout
    elif 'w' in mode:
        return subprocess.Popen(command, shell=True, bufsize=buff,
                                stdin=subprocess.PIPE).stdin
    return None

NORMAL = 0
PROCESS = 1
PARALLEL = 2
WHICH_GZIP = which("gzip")
WHICH_PIGZ = which("pigz")

def open_gz(filename, mode='r', buff=1024*1024, external=PARALLEL):
    if external == None or external == NORMAL:
        import gzip
        return gzip.GzipFile(filename, mode, buff)
    elif external == PROCESS:
        if not WHICH_GZIP:
            return open_gz(filename, mode, buff, NORMAL)
        if 'r' in mode:
            return open_pipe("gzip -dc " + filename, mode, buff)
        elif 'w' in mode:
            return open_pipe("gzip >" + filename, mode, buff)
    elif external == PARALLEL:
        if not WHICH_PIGZ:
            return open_gz(filename, mode, buff, PROCESS)
        if 'r' in mode:
            return open_pipe("pigz -dc " + filename, mode, buff)
        elif 'w' in mode:
            return open_pipe("pigz >" + filename, mode, buff)
    return None

def zopen(filename, mode='r', buff=1024*1024, external=PARALLEL):
    """
    Open pipe, zipped, or unzipped file automagically
    # external == 0: normal zip libraries
    # external == 1: (zcat, gzip) or (bzcat, bzip2)
    # external == 2: (pigz -dc, pigz) or (pbzip2 -dc, pbzip2)
    """
    if 'r' in mode and 'w' in mode:
        return None
    if filename.startswith('!'):
        return open_pipe(filename[1:], mode, buff)
    elif filename.endswith('.gz'):
        return open_gz(filename, mode, buff, external)
    else:
        return open(filename, mode, buff)
    return None

