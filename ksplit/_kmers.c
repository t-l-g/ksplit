#include <inttypes.h>

const uint8_t reverse4table[] = {
    0,  8,  4, 12,  2, 10,  6, 14,  1,  9,  5, 13,  3, 11,  7, 15
};

uint8_t reverse8(uint8_t x) {
    return (reverse4table[x & 0xf] << 4)
           | reverse4table[x >> 4];
}

uint32_t reverse32(uint32_t v) {
    return (reverse8(v & 0xff)  << 24)
        | (reverse8((v >> 8) & 0xff) << 16)
        | (reverse8((v >> 16) & 0xff) << 8)
        | reverse8((v >> 24) & 0xff);
}
uint64_t reverse64(uint64_t v) {
    return ((uint64_t)reverse32(v & 0xffffffff)  << 32)
        | reverse32(v >> 32);
}

uint64_t canonical_kmer(uint64_t kmer) {
    uint64_t alternative = reverse64(kmer) >> 2;
    if (alternative < kmer) return alternative;
    return kmer;
}

uint64_t encode_nt(char nt) {
    if (nt == 'A') return 0;
    if (nt == 'C') return 1;
    if (nt == 'T') return 2;
    if (nt == 'G') return 3;
    return -1;

}
