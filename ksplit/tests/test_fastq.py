from ksplit.fastq import fastq_iter, _compatible
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
    assert x.seq == 'TATATATATATTTCTTGTAATTTGTTGGAATACGAGAACATCGTCAATAATATATCGTATGAATTGAACCACACGGCACA'
    [[x,x2]] = list(fastq_iter(StringIO(fq+fq)))
    assert x.seq == 'TATATATATATTTCTTGTAATTTGTTGGAATACGAGAACATCGTCAATAATATATCGTATGAATTGAACCACACGGCACA'
    assert x2.seq == 'TATATATATATTTCTTGTAATTTGTTGGAATACGAGAACATCGTCAATAATATATCGTATGAATTGAACCACACGGCACA'

def test_interleaved_fastq():
    with open(path.join(BASEDIR, 'interleaved.fq'), 'r') as ifile:
        lens = [len(e)
                for e in fastq_iter(ifile)]
    assert lens == [2,2,1,1]

def test_compatible():
    assert _compatible('@SRR4052021.7 7/1', '@SRR4052021.7 7/2')
    assert not _compatible('@SRR4052021.7 7/1', '@SRR4052021.7 8/1')
    assert _compatible('@SRR4052021.7', '@SRR4052021.7')
    assert not _compatible('@SRR4052021.7 xz', '@SRR4052021.7 xy')
