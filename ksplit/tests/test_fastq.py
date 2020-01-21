from ksplit.fastq import fastq_iter
from os import path

BASEDIR = path.join(path.dirname(path.abspath(__file__)), 'data')

fq = '''\
@ReadP
TATATATATATTTCTTGTAATTTGTTGGAATACGAGAACATCGTCAATAATATATCGTATGAATTGAACCACACGGCACA
+
HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
'''

def test_iter_fastq():
    from io import StringIO
    [[x]] = list(fastq_iter(StringIO(fq)))
    assert x == 'TATATATATATTTCTTGTAATTTGTTGGAATACGAGAACATCGTCAATAATATATCGTATGAATTGAACCACACGGCACA'
    [[x,x2]] = list(fastq_iter(StringIO(fq+fq)))
    assert x == 'TATATATATATTTCTTGTAATTTGTTGGAATACGAGAACATCGTCAATAATATATCGTATGAATTGAACCACACGGCACA'
    assert x2 == 'TATATATATATTTCTTGTAATTTGTTGGAATACGAGAACATCGTCAATAATATATCGTATGAATTGAACCACACGGCACA'

def test_interleaved_fastq():
    with open(path.join(BASEDIR, 'interleaved.fq'), 'r') as ifile:
        lens = [len(e)
                for e in fastq_iter(ifile)]
    assert lens == [2,2,1,1]
