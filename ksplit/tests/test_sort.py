import numpy as np
from io import BytesIO
from ksplit import ksplit, sort
import numpy as np
from os import stat

fq = '''@ReadP
TATATATATATTTCTTGTAATTTGTTGGAATACGAGAACATCGTCAATAATATATCGTATGAATTGAACCACACGGCACA
+
HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
'''

def test_sort_partials(tmpdir):
    import operator
    import functools

    out = BytesIO()
    ksplit.encode_fastq({}, BytesIO(fq.encode('ascii')), out)
    data = out.getvalue()

    full = sort._from_buffer(data)
    sp = sort.sort_partials({}, BytesIO(data), tmpdir, block_nbytes=16*8*2)
    assert sum(stat(f).st_size for f in sp) == len(data)

    blocks = [sort._from_buffer(open(f, 'rb').read()) for f in sp]
    assert functools.reduce(operator.or_, (set(block.T[0]) for block in blocks)) == set(full.T[0])

    for b in blocks:
        p = b.T[0].copy()
        p.sort()
        assert np.all(p == b.T[0])
