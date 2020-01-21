import pyximport
import numpy as np
pyximport.install(setup_args={
    'include_dirs': np.get_include()
    })

from . import kmers
from .fastq import fastq_iter


def encode_fastq(args, ifile, out):
    '''
    Returns the number of sequences encoded

    Example
    -------
    encode_fastq( { .. }, open(ifname, 'rb'), open(ofname, 'wb'))
    '''
    for i,seq in enumerate(fastq_iter(ifile)):
        seq = seq.replace(b'N', b'A')
        ks = kmers.kmers(seq)
        if ks is None:
            raise ValueError("Something wrong!")
        ixs = np.repeat(np.array([i], dtype=np.uint64), len(ks))
        encoded = np.vstack([ks, ixs]).T.ravel()
        out.write(encoded.data)
    return i + 1
