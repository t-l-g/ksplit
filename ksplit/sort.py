from os import path, makedirs
import numpy as np

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
