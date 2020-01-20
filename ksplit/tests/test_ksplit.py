from io import BytesIO
import ksplit.ksplit
import numpy as np

fq = '''@ReadP
TATATATATATTTCTTGTAATTTGTTGGAATACGAGAACATCGTCAATAATATATCGTATGAATTGAACCACACGGCACA
+
HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
'''

def test_encode_fastq():
    out = BytesIO()
    ksplit.ksplit.encode_fastq({}, BytesIO(fq.encode('ascii')), out)
    assert set(np.frombuffer(out.getvalue(), np.uint64).reshape((-1, 2)).T[1]) == set([0])

    out = BytesIO()
    ksplit.ksplit.encode_fastq({}, BytesIO((fq+fq).encode('ascii')), out)
    assert set(np.frombuffer(out.getvalue(), np.uint64).reshape((-1, 2)).T[1]) == set([0, 1])
