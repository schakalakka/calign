import argparse
import glob
import os
import time

from Bio import SeqIO
from Bio.Alphabet import generic_dna
from cffi import FFI

BASE_DIR = '/home/andreas/GDrive/workspace/sparsedata/'
FILES = sorted(glob.glob(BASE_DIR + '*/*.fasta'))

ffi = FFI()
lib = ffi.dlopen('./lib/libcalignshared.so')
print('Loaded lib {}'.format(lib))
ffi.cdef('''
void alignment_score_all(const char * filename, int len, const char** strings, char* alignment_type);
''')
ffi.cdef('''
void banded_alignment_score_all(const char * filename, int len, const char** strings, int lower_diag, int upper_diag);
''')


# ffi.cdef('''typedef struct {
#     char *a;
#     unsigned int alen;
#     char *b;
#     unsigned int blen;
# } seq_pair;
# typedef seq_pair *seq_pair_t;''')
#
#
#
# ffi.cdef('''int alignment_score(seq_pair_t problem, char* alignment_type);''')
# ffi.cdef('''int semiglobal_alignment_score(seq_pair_t problem, int, int, int, int);''')
# ffi.cdef('''int alignment_score2(const char* , int , const char* , int );''')


def run(dir, ldiag, udiag):
    fasta_file = glob.glob(dir + '*.fasta')[0]
    print(fasta_file)
    reads = [ffi.new('const char[]', (bytes(str(seqRecord.seq), 'ascii'))) for seqRecord in
           SeqIO.parse(fasta_file, "fasta", generic_dna)]
    t0 = time.time()
    sub_dir = dir.split('/')[-2]
    if ldiag == udiag == None:
        filename = ffi.new('const char[]', bytes(sub_dir + 'calign.score', 'ascii'))
        lib.alignment_score_all(filename, len(reads), reads, b'semiglobal')
        t1 = time.time()
        os.system('mv {} {}'.format(sub_dir + 'calign.score', dir + 'calign.score'))
        with open(dir + 'calign.time', 'w') as f:
            f.write('{}'.format(t1 - t0))
    else:
        filename = ffi.new('const char[]', bytes(sub_dir + 'calign_{}_{}.score'.format(ldiag, udiag), 'ascii'))
        lib.banded_alignment_score_all(filename, len(reads), reads, ldiag, udiag)
        t1 = time.time()
        os.system('mv {} {}'.format(sub_dir + 'calign_{}_{}.score'.format(ldiag, udiag),
                                    dir + 'calign_{}_{}.score'.format(ldiag, udiag)))
        with open(dir + 'calign_{}_{}.time'.format(ldiag, udiag), 'w') as f:
            f.write('{}'.format(t1 - t0))


# def get_score(a_str, b_str):
#     a = bytes(a_str, 'ascii')
#     b = bytes(b_str, 'ascii')
#     lib.alignment_score2(ffi.new('const char[]', a), len(a), ffi.new('const char[]', b), len(b))
#     foo = ffi.new('seq_pair_t')
#     foo.a = ffi.new('const char[]', a)
#     foo.alen = len(a)
#     foo.b = ffi.new('const char[]', b)
#     foo.blen = len(b)
#     bla = ffi.new('const char[]', b'semiglobal')
#     bar = lib.alignment_score(foo, bla)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('dir', help='Directory of reads library file', type=str)
    parser.add_argument('-l', '--lower', help='Lower diagonal', type=int, default=None)
    parser.add_argument('-u', '--upper', help='Upper diagonal', type=int, default=None)
    # parser.add_argument('-m', '--match', help='Match score', type=int, default=1)
    # parser.add_argument('-n', '--mismatch', help='Mismatch penalty', type=int, default=-5000)
    # parser.add_argument('-g', '--gap', help='Gap penalty', type=int, default=-5000)

    args = parser.parse_args()
    run(args.dir, args.lower, args.upper)

    # if len(sys.argv) == 2:
    #     run(sys.argv[1])
    # else:
    #     print('Usage: calignrun.py FASTA_DIR')
    #     sys.exit(1)
