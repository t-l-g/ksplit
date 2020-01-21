from os import path, makedirs
import numpy as np

class file_buffer:
    '''file object wrapper with cheap peek()ing'''
    def __init__(self, fname):
        self.handle = open(fname, 'rb')
        self.buf = b''


    def close(self):
        self.handle.close()


    def __del__(self):
        self.close()


    def peek(self, nbytes):
        if nbytes > len(self.buf):
            nbytes -= len(self.buf)
            read_size = max(4096, nbytes)
            self.buf += self.handle.read(read_size)
        return self.buf[:nbytes]


    def consume(self, nbytes):
        self.buf = self.buf[nbytes:]


def _from_buffer(data):
    return np.frombuffer(data, np.uint64).reshape((-1, 2))

def _read_blocks(ifile, nbytes):
    n = nbytes / (2 * 8)
    while True:
        data = ifile.read(nbytes)
        if not data:
            return
        yield _from_buffer(data).copy()

def sort_partials(args, encoded_stream, tdir, *, block_nbytes=1024*1024*1024):
    splits_dir = path.join(tdir, 'splits')
    makedirs(splits_dir)
    partials = []
    for block in _read_blocks(encoded_stream, block_nbytes):
        ofname = path.join(splits_dir, 'split_{}'.format(len(partials)))
        partials.append(ofname)
        with open(ofname, 'wb') as out:
            block.sort(axis=0)
            out.write(block.data)
    return partials


def merge_stream(bufs, out, block_nbytes):
    while bufs:
        min_k = np.uint64(-1)
        # For the first pass, we could, in principle, just seek() to the right
        # position and read those 8 Bytes
        for b in bufs:
            ch = _from_buffer(b.peek(block_nbytes))
            if not ch.size:
                continue
            cur = ch.T[0][-1]
            if cur < min_k: min_k = cur
        nbufs = []
        cur = []
        for b in bufs:
            ch = _from_buffer(b.peek(block_nbytes))
            if not ch.size:
                continue
            cur_c = np.searchsorted(ch.T[0], min_k, 'right')
            b.consume(cur_c * 8 * 2)
            cur.append(ch[:cur_c])
            nbufs.append(b)
        bufs = nbufs
        if cur:
            cur = np.concatenate(cur)
            cur.sort(axis=0)
            out.write(cur.data)
