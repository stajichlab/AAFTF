# this module runs contamination screening

import sys, os
import shutil

from subprocess import call, run, Popen, PIPE, STDOUT

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna


