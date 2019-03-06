#ifndef BFMEM_H
#define BFMEM_H

#include <stdint.h>
#include <stdio.h>
#include <sys/types.h>
#include "khash.h"
#include "kvec.h"
#include "ksort.h"

struct mm128_t{
	uint64_t x, y;
	mm128_t(): x(0), y(0) {}
	mm128_t(uint64_t _x, uint64_t _y): x(_x), y(_y) {}
};

typedef struct { size_t n, m; mm128_t *a; } mm128_v;

typedef struct {
	mm128_v a;   // (hash value, position) array
	int32_t n;   // size of the _p_ array
	uint64_t *p; // position array for hash value appearing >1 times
	void *h;     // hash table indexing _p_ and hash value appearing once
} mm_idx_bucket_t;

typedef struct {
	uint32_t n, b;  // number of reference sequences
	mm_idx_bucket_t *B;
} mm_idx_t;

mm_idx_t *mm_idx_init(int b) {
	mm_idx_t *mi;
	mi = (mm_idx_t*)calloc(1, sizeof(mm_idx_t));
	mi->B = (mm_idx_bucket_t*)calloc(1<<b, sizeof(mm_idx_bucket_t));
	mi->b = b;
	return mi;
}

void radix_sort_128x(mm128_t *beg, mm128_t *end);

/*********************************
 * Sort and generate hash tables *
 *********************************/

#define idx_hash(a) ((a)>>1)
#define idx_eq(a, b) ((a)>>1 == (b)>>1)
KHASH_INIT(idx, uint64_t, uint64_t, 1, idx_hash, idx_eq)
typedef khash_t(idx) idxhash_t;

/* this function is https://github.com/lh3/minimap/blob/master/index.c
 * authors: Heng Li
 * changes: no
 */
const uint64_t *mm_idx_get(const mm_idx_t *mi, uint64_t minier, int *n) {
	int mask = (1<<mi->b) - 1;
	khint_t k;
	mm_idx_bucket_t *b = &mi->B[minier&mask];
	idxhash_t *h = (idxhash_t*)b->h;
	*n = 0;
	if (h == 0) return 0;
	k = kh_get(idx, h, minier>>mi->b<<1);
	if (k == kh_end(h)) return 0;
	if (kh_key(h, k)&1) {
		*n = 1;
		return &kh_val(h, k);
	} else {
		*n = (uint32_t)kh_val(h, k);
		return &b->p[kh_val(h, k)>>32];
	}
}

/* this function is https://github.com/lh3/minimap/blob/master/index.c
 * authors: Heng Li
 * changes: delete some lines
 */
void mm_idx_destroy(mm_idx_t *mi) {
	int i;
	if (mi == 0) return;
	for (i = 0; i < (1<<mi->b); ++i) {
		free(mi->B[i].p);
		free(mi->B[i].a.a);
		kh_destroy(idx, (idxhash_t*)mi->B[i].h);
	}
	free(mi->B);
	free(mi);
}

// 2 4 8 16 32 64 128 at most 128 files
// inline int hash2(const size_t &a, const size_t &b, const size_t &c) {
inline int hash2(const size_t &b) {
	// return (a >> 7)&1;
	return (b >> 9)&1;
}

// inline int hash4(const size_t &a, const size_t &b, const size_t &c) {
	// return (((a >> 7)&1)<<1) + ((b >> 5)&1);
inline int hash4(const size_t &b) {
	return (((b >> 9)&1)<<1) + ((b >> 5)&1);
}

// inline int hash8(const size_t &a, const size_t &b, const size_t &c) {
inline int hash8(const size_t &b) {
	return (((b >> 9)&1)<<2) + (((b >> 5)&1)<<1) + ((b >> 2)&1);
}

// inline int hash16(const size_t &a, const size_t &b, const size_t &c) { //2^4; 2 1 1
inline int hash16(const size_t &b) {
	return (((b >> 9)&3)<<2) + (((b >> 5)&1)<<1) + ((b >> 2)&1);
}

// inline int hash32(const size_t &a, const size_t &b, const size_t &c) { //2^5; 2 2 1
inline int hash32(const size_t &b) {
	return (((b >> 9)&3)<<3) + (((b >> 5)&3)<<1) + ((b >> 2)&1);
}

// inline int hash64(const size_t &a, const size_t &b, const size_t &c) { //2^6; 3 2 1
inline int hash64(const size_t &b) {
	return (((b >> 9)&7)<<3) + (((b >> 5)&3)<<1) + ((b >> 2)&1);
}

// inline int hash128(const size_t &a, const size_t &b, const size_t &c) { //2^7; 3 3 1
inline int hash128(const size_t &b) {
	return (((b >> 12)&7)<<4) + (((b >> 5)&7)<<1) + ((b>>2)&1);
}


#endif