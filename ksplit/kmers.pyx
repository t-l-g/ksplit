import numpy as np
cimport numpy as np

cdef extern from "_kmers.c":
    long long encode_nt(char);
    long long canonical_kmer(long long);

cdef _kmers(char* seq, int n):
    cdef int KMER_SIZE = 31
    out = np.zeros(n - KMER_SIZE + 1, np.uint64)

    cdef long long kmer = 0
    cdef long long nte = 0
    cdef int j = 0
    for i in range(n):
        kmer >>= 2
        nte = encode_nt(seq[i])
        if nte == -1:
            return
        kmer |= nte << ((KMER_SIZE - 1) * 2);
        if i >= (KMER_SIZE - 1):
            out[j] = canonical_kmer(kmer)
            j += 1
    return out

def kmers(seq):
    return _kmers(<bytes>seq, len(seq))
