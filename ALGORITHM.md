# K-split algorithm

## High level

1. Compute kmer profiles & sort
2. Run union-find to find connected components
3. If necessary, break up large CCs using min-cut
4. Split the resulting breaks into blocks

All steps must be done in limited memory (out-of-core)

## Kmer profile computation and sorting

## Union-find

## Min-cut

The [...] algorithm is heuristic, but very easy to implement in an out-of-core fashion:

1. Start with V components
2. Randomly sort all the edges in the graph
3. Go through the (randomly sorted) edge list and merge them until you have two
   groups. With high expectation, you will have a good min-cut approximation




