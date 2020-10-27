# cython: language_level=3
import numpy as np
cimport numpy as np

from libc.stdint cimport uint8_t, uint64_t

cdef uint64_t encode_nt(char nt) nogil:
    if nt == b'A': return 0
    if nt == b'C': return 1
    if nt == b'T': return 2
    if nt == b'G': return 3
    return -1

cdef uint64_t encode_nt_c(char nt) nogil:
    if nt == b'A': return encode_nt(b'T')
    if nt == b'C': return encode_nt(b'G')
    if nt == b'T': return encode_nt(b'A')
    if nt == b'G': return encode_nt(b'C')
    return -1

cdef _kmers(char* seq, int n):
    cdef int KMER_SIZE = 31
    out = np.zeros(n - KMER_SIZE + 1, np.uint64)

    cdef uint64_t kmer = 0
    cdef uint64_t kmer_rc = 0
    cdef uint64_t nte = 0
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
    '''Compute all kmers for input nucleotide sequence

    Parameters
    ----------
    seq : bytes
        input nucleotide sequence

    Returns
    -------
    kmers : ndarray
        kmers encoded as integers in a NumPy array
    '''
    return _kmers(<bytes>seq, len(seq))

