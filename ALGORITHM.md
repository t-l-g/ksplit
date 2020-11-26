# K-split algorithm

## High level

1. Compute kmer profiles & sort
2. Run union-find to find connected components
3. If necessary, break up large CCs using min-cut
4. Split the resulting breaks into blocks

All steps must be done in limited memory (out-of-core)

## Kmer profile computation and sorting

### On disk format

For sorting, the kmers are saved on disk. They are saved in a binary format
such that each consecutive 8-Bytes (64 bits) correspond to an integer (machine
format is assumed as these will not be shared across a network) and each
consecutive two numbers form a pair (_i.e._, `[4302, 0, 9721, 0, 231, 0]`
should be interpreted as `[(4302, 0), (9721, 0), (231, 0)]`). The first element
of the pair is the encoded kmer, while the second element is the sequence
identifier (sequences being identified by integer IDs, starting at zero).

## Union-find

## Min-cut

The [...] algorithm is heuristic, but very easy to implement in an out-of-core fashion:

1. Start with V components
2. Randomly sort all the edges in the graph
3. Go through the (randomly sorted) edge list and merge them until you have two
   groups. With high expectation, you will have a good min-cut approximation




