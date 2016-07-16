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


def run(file):
    foo = [ffi.new('const char[]', (bytes(str(seqRecord.seq), 'ascii'))) for seqRecord in
           SeqIO.parse(file, "fasta", generic_dna)]
    t0 = time.time()
    lib.alignment_score_all(len(foo), foo, b'semiglobal')
    t1 = time.time()
    os.system('mv out.out {}'.format(file.split('.fasta')[0] + '.calignscore'))
    with open(file.split('.fasta')[0] + '.cact', 'w') as f:
        f.write('{}'.format(t1 - t0))


if len(sys.argv) == 1:
    for file in FILES:
        print(file.split('/')[-2])
        run(file)
else:
    foo = glob.glob(sys.argv[1] + '*.fasta')
    run(glob.glob(sys.argv[1] + '*.fasta')[0])
