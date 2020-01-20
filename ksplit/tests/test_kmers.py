import pyximport
pyximport.install(build_dir='.', build_in_temp=False)


from ksplit import kmers

def test_kmers():
    testing = 'ATTTAACATGAGATAACATGCATGCATGCATTGCGGCTCAGCTAGTCAGCTAGCTAGCTAGCTACGATCGATCGTAGCATCGATCGATCGATCGATCGATCGATCGTACGTACGTAGCTACGATCGTAGCTAGCTAG'

    testing = testing.encode('ascii')
    print(kmers.kmers(testing))

    testing = 'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'

    testing = testing.encode('ascii')
    print(kmers.kmers(testing))


def test_error_detection():
    testing = 'ATTTAACATGAGXTAACATGCATGCATGCAT'
    testing = testing.encode('ascii')
    assert kmers.kmers(testing) is None


def test_kmers1():
    testing = 'TTTCTTTTTTTTTTTTTTTTTTTTTTTTTTT'
    testing = testing.encode('ascii')
    ks = kmers.kmers(testing)
    assert len(testing) == 31
    assert len(ks) == 1

