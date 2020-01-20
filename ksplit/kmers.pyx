import numpy as np
cimport numpy as np

cdef long long encode_nt(char nt):
    if nt == 'A': return 0
    if nt == 'C': return 1
    if nt == 'T': return 2
    if nt == 'G': return 3
    return -1

cdef long long encode_nt_c(char nt):
    if nt == 'A': return encode_nt('T')
    if nt == 'C': return encode_nt('G')
    if nt == 'T': return encode_nt('A')
    if nt == 'G': return encode_nt('C')
    return -1

cdef _kmers(char* seq, int n):
    cdef int KMER_SIZE = 31
    out = np.zeros(n - KMER_SIZE + 1, np.uint64)

    cdef long long kmer = 0
    cdef long long kmer_rc = 0
    cdef long long nte = 0
    cdef int j = 0
    for i in range(n):
        # The kmer of the reverse complement should be the same as that of the
        # string. So, we compute both and always return the minimum

        kmer >>= 2
        kmer_rc <<= 2
        kmer_rc &= ~(0x3 << (2*(KMER_SIZE-1)))
        nte = encode_nt(seq[i])
        if nte == -1:
            return
        kmer |= nte << ((KMER_SIZE - 1) * 2);
        kmer_rc |= encode_nt_c(seq[i])
        if i >= (KMER_SIZE - 1):
            out[j] = min(kmer, kmer_rc)
            j += 1
    return out

def kmers(seq):
    return _kmers(<bytes>seq, len(seq))
