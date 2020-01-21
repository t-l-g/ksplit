import numpy as np
cimport numpy as np

from libc.stdint cimport uint32_t, uint64_t

cdef uint32_t _follow(uint32_t [::1] uf, uint32_t v) nogil:
    while uf[v] != v:
        v = uf[v]
    return v

cdef void _link(uint32_t [::1] uf, uint32_t a, uint32_t b) nogil:
    cdef uint32_t c = _follow(uf, b)
    uf[b] = c
    uf[a] = c

cdef void _make_canonical(uint32_t [::1] uf, int n) nogil:
    for i in range(n):
        uf[i] = _follow(uf, i)

def union_find(rs, int n, int block_nbytes):
    cdef uint64_t pkm = -1
    cdef uint32_t pn = 0
    cdef uint64_t km
    cdef uint32_t r
    cdef uint64_t [:, ::1] ch_view
    cdef int i, ch_len

    uf_arr = np.arange(n, dtype=np.uint32)
    cdef uint32_t [::1] uf = uf_arr

    while True:
        data = rs.read(block_nbytes)
        if not data:
            break
        ch = np.frombuffer(data, np.uint64).reshape((-1, 2)).copy()
        ch_len = len(ch)
        ch_view = ch
        for i in range(ch_len):
            km = ch_view[i, 0]
            r = ch_view[i, 1]
            if km == pkm:
                _link(uf, r, pn)
            else:
                pkm = km
                pn = r
    _make_canonical(uf, n)

    return uf_arr
