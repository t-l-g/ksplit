import pyximport
import numpy as np
pyximport.install(setup_args={
    'include_dirs': np.get_include()
    })

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

def rc(t):
    rcd = {
        'A': 'T',
        'T': 'A',
        'C': 'G',
        'G': 'C',
    }
    return ''.join([rcd[c] for c in t[::-1]])


def test_kmers_reverse():
    for t in [
            'TTATACATACTGTTGGTATGATAATAGTATA',
            'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT',
            'ATTTTTTTTTTTTTTTTTTTTTTTTTTTTTT',
            'AAAAAAAAAAATTTTTTTTTTTTTTTTTTTT',
            ]:
        assert np.all(kmers.kmers(t.encode('ascii'))
                    == kmers.kmers(rc(t).encode('ascii')))


def test_kmers_reverse_embed():
    k = 'TTATACATACTGTTGGTATGATAATAGTATA'

    t0 = k + 'C'
    t1 = rc(k) + 'T'
    assert kmers.kmers(t0.encode('ascii'))[0] == kmers.kmers(t1.encode('ascii'))[0]

def test_max62():
    assert len('{:b}'.format(kmers.kmers('GACATAGCGACGCGGACCCCCTTTTTTTTTTGG'.encode('ascii')).max())) <= 62
