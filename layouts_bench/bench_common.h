/* =========================================================================
 * bench_common.h
 *
 * Shared types, timing utilities, Peano-Hilbert key generation, radix sort,
 * permutation helpers, and PRNG for the OG3 data-layout benchmarks.
 * ========================================================================= */

#ifndef BENCH_COMMON_H
#define BENCH_COMMON_H

#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif
#ifndef _ISOC11_SOURCE
#define _ISOC11_SOURCE
#endif

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#ifndef OG3_ATOM_ALIGN
#define OG3_ATOM_ALIGN 64
#endif

#define OG3_ATOM __attribute__((aligned(OG3_ATOM_ALIGN)))

/* Scalar aliases. */
typedef double   pos_t;
typedef double   acc_t;
typedef float    real_t;
typedef uint64_t pkey_t;
typedef uint64_t ppid_t;
typedef int64_t  time_int_t;

/* Particle counts and particle-array indices. */
typedef uint64_t count_t;

/* Particle types. */
typedef enum {
    PT_GAS  = 0,
    PT_DM1  = 1,
    PT_DM2  = 2,
    PT_DM3  = 3,
    PT_STAR = 4,
    PT_BH   = 5,
    NTYPES  = 6
} part_type_t;

#define PT_IS_DM(t) ((t) == PT_DM1 || (t) == PT_DM2 || (t) == PT_DM3)

#ifndef NMET
#define NMET 16
#endif

typedef struct {
    count_t n_gas;
    count_t n_dm;
    count_t n_star;
    count_t n_bh;
} pcount_t;

void counts_from_total(count_t N_total, pcount_t * restrict out);

static inline count_t total_parts(const pcount_t * restrict c)
{
    return c->n_gas + c->n_dm + c->n_star + c->n_bh;
}

/* Per-process CPU timer. */
static inline double cpu_time_p(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

#define TIME_BLOCK(label, stmt)                                             \
    do {                                                                    \
        double _t0 = cpu_time_p();                                          \
        stmt;                                                               \
        double _t1 = cpu_time_p();                                          \
        fprintf(stderr, "[time] %-32s %10.3f ms\n",                       \
                (label), 1000.0 * (_t1 - _t0));                             \
    } while (0)

#define OG3_ALLOC_ALIGN 64

static inline void *og3_aligned_alloc(size_t bytes)
{
    if (bytes == 0) return NULL;
    size_t rounded = (bytes + (OG3_ALLOC_ALIGN - 1)) &
                     ~((size_t)OG3_ALLOC_ALIGN - 1);
    return aligned_alloc(OG3_ALLOC_ALIGN, rounded);
}

/* Peano-Hilbert keys. ph_key_from_ijk expects x,y,z < 2^21. */
pkey_t ph_key_from_ijk(uint32_t i, uint32_t j, uint32_t k);
void gen_keys(pkey_t * restrict keys, count_t N,
              uint32_t n_cells, uint64_t seed);

void assign_types(uint8_t * restrict types,
                  const pcount_t * restrict c, uint64_t seed);

/* Stable LSD radix sort by 64-bit key. Payload is count_t so permutations
 * remain valid beyond 2^32 particles per rank. */
void radix_sort_u64_count(const pkey_t  * restrict keys_in,
                          const count_t * restrict pl_in,
                          pkey_t        * restrict keys_out,
                          count_t       * restrict pl_out,
                          count_t N);

void invert_permutation(const count_t * restrict perm,
                        count_t       * restrict inv,
                        count_t N);

int verify_sort(const pkey_t * restrict keys_sorted, count_t N);

void gen_sparse_permutation(count_t * restrict perm, count_t N,
                            double fraction, uint64_t seed);

typedef struct {
    count_t i, j;
} swap_pair_t;

void gen_swap_pairs_same_type(swap_pair_t   * restrict pairs,
                              count_t n_pairs,
                              const uint8_t * restrict types,
                              count_t N,
                              uint64_t seed);

count_t count_permutation_changed(const count_t * restrict perm, count_t N);

/* Apply a gather permutation: after the call, array[i] contains the old
 * array[perm[i]]. The helper chooses a dense gather-copy-back path for dense
 * permutations and a changed-destination-only path for sparse permutations.
 * tmp_bytes must hold N * elem_size bytes. changed_idx must hold N entries. */
void og3_apply_perm_bytes(void * restrict array,
                          const count_t * restrict perm,
                          count_t N,
                          size_t elem_size,
                          void * restrict tmp_bytes,
                          count_t * restrict changed_idx);

size_t og3_perm_dense_copy_bytes(count_t N, size_t elem_size);
size_t og3_perm_sparse_copy_bytes(count_t changed, size_t elem_size);

#ifdef DEBUG
uint64_t og3_debug_marker64(ppid_t id, uint64_t tag);
uint32_t og3_debug_marker32(ppid_t id, uint32_t tag);
#endif

/* PRNG. */
typedef struct { uint64_t s[4]; } rng_t;

void     rng_seed(rng_t * restrict r, uint64_t seed);
uint64_t rng_next(rng_t * restrict r);
double   rng_uniform(rng_t * restrict r);

#endif /* BENCH_COMMON_H */
