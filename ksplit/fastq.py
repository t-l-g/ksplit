def fastq_iter(ifile):
    '''
    for seq in fastq_iter(open('example.fq', 'rb'))):
        print(seq)
    '''
    while True:
        h = ifile.readline()
        seq = ifile.readline()
        _ = ifile.readline()
        qs = ifile.readline()
        if not qs:
            break
        yield seq.strip()

