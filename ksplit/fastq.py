def fastq_iter1(ifile):
    '''
    for h, seq in fastq_iter1(open('example.fq', 'rb'))):
        print(h, seq)
    '''
    while True:
        h = ifile.readline()
        seq = ifile.readline()
        _ = ifile.readline()
        qs = ifile.readline()
        if not qs:
            break
        yield h, seq.strip()


def fastq_iter(base):
    '''
    Read in an interleaved FastQ

    This is an iterator over **just the sequence information** and supports
    optionally-interleaved format. It yields tuples of 1 or 2 sequences.

    Two adjacent reads are considered to form a pair if their headers are
    identical; otherwise they are considered single-end.

    Example
    -------

    for seqs in fastq_iter(open('interleaved.fq', 'rb')):
        if len(seqs) == 1:
            [s] = seqs
        else:
            [mate1, mate2] = seqs
    '''

    prev_h = ''
    prev_seq = None
    for h,seq in fastq_iter1(base):
        if prev_h != h:
            if prev_seq is not None:
                yield (prev_seq,)
            prev_h = h
            prev_seq = seq
        else:
            yield (prev_seq, seq)
            prev_h = ''
            prev_seq = None
    if prev_seq is not None:
        yield (prev_seq,)

