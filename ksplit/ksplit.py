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


    Check the file `ALGORITHM.md` for a description of the on-disk format

    Parameters
    ----------

    args : arguments
    ifile (file-like): input file object (should be a FastQ file)
    ofile (file-like): output file object (should be opened in binary mode)

    Example
    -------
    encode_fastq( { .. }, open(ifname, 'rb'), open(ofname, 'wb'))
    '''
    for i,seqs in enumerate(fastq_iter(ifile)):
        ks = []
        for seq in seqs:
            # Below is a major hack, but this was actually done by SOAP too
            seq = seq.seq.replace(b'N', b'A')
            k = kmers.kmers(seq)
            if k is None:
                raise ValueError("Something wrong!")
            ks.append(k)

        # TODO This reorganization of the ks could be done in one step without
        # generating the intermediate arrays, but it would require somewhat
        # more complex code:
        if len(ks) == 1:
            [ks] = ks
        else:
            ks = np.concatenate(ks)
        ixs = np.repeat(np.array([i], dtype=np.uint64), len(ks))
        encoded = np.vstack([ks, ixs]).T.ravel()
        out.write(encoded.data)
    return i + 1
