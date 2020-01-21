import pyximport
import numpy as np
pyximport.install(setup_args={
    'include_dirs': np.get_include()
    })

from ksplit import graph

def test_union_find():
    from io import BytesIO

    # The graph below should have 3 CCs:
    # 0-1-3-6-7
    # 2-5
    # 4

    r = np.array([
       [0,0],
       [0,1],

       [1,1],
       [1,3],
       [1,6],
       [1,7],

       [2,4],

       [3,0],
       [3,6],

       [4,2],

       [5,5],
       [5,2],

       [6,4],

       [7,0],

       [8,1]], dtype=np.uint64)

    uf = graph.union_find(BytesIO(r.data), 8, 32)

    assert len(set(uf[[0,1,6,7]]))   == 1
    assert len(set(uf[[2,5]]))       == 1
    assert len(set(uf[[2,5, 0]]))    == 2
    assert len(set(uf[[2,5, 0, 4]])) == 3
    assert len(set(uf)) == 3
