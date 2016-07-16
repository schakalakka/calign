import glob
import os
import sys
import time

from Bio import SeqIO
from Bio.Alphabet import generic_dna
from cffi import FFI

BASE_DIR = '/home/andreas/GDrive/workspace/sparsedata/'
FILES = sorted(glob.glob(BASE_DIR + '*/*.fasta'))

ffi = FFI()
lib = ffi.dlopen('./lib/libsharedcalign.so')
print('Loaded lib {}'.format(lib))
ffi.cdef('''
void alignment_score_all(int len, const char** strings, char* alignment_type);
''')


def run(DIR):
    fasta_file = glob.glob(DIR + '*.fasta')[0]
    foo = [ffi.new('const char[]', (bytes(str(seqRecord.seq), 'ascii'))) for seqRecord in
           SeqIO.parse(fasta_file, "fasta", generic_dna)]
    t0 = time.time()
    lib.alignment_score_all(len(foo), foo, b'semiglobal')
    t1 = time.time()
    os.system('mv out.out {}'.format(DIR + 'calign.score'))
    with open(DIR + 'calign.time', 'w') as f:
        f.write('{}'.format(t1 - t0))


if len(sys.argv) == 2:
    run(sys.argv[1])
else:
    print('Usage: calignrun.py FASTA_DIR')
    sys.exit(1)
