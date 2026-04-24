#ifndef BOXLEAVES_H
#define BOXLEAVES_H

#include "bench_common.h"

#define BOXLEAF_INVALID ((count_t)~(count_t)0)
#define OG3_PH_BITS 21u
#define OG3_PH_DIM 3u
#define OG3_PH_KEY_BITS (OG3_PH_BITS * OG3_PH_DIM)

typedef struct {
    /* PH cube represented by [key_begin, key_end). */
    pkey_t key_begin;
    pkey_t key_end;

    /* Current contiguous common-array range occupied by this leaf. */
    count_t begin;
    count_t end;

    /* Current gas sub-range inside [begin,end). The partition benchmark keeps
     * gas first within each leaf, so [gas_begin,gas_end) is contiguous. */
    count_t gas_begin;
    count_t gas_end;

    /* User-code-compatible description of the PH cube. octant is key_begin,
     * masked at level. level uses the natural 0..21 convention here: level 0
     * is the full box, level 21 is one integer PH cell per axis. */
    pkey_t octant;
    uint32_t side_cells;
    uint8_t level;
} boxleaf_t;

typedef struct {
    boxleaf_t *leaves;
    count_t n_leaf;
    count_t capacity;

    count_t target;
    count_t min_particles;
    count_t max_particles;
    double tolerance;
    uint32_t max_side_cells;

    count_t underfull;
    count_t overfull;
    count_t spatially_forced;
} boxleaf_set_t;

void boxleaf_set_free(boxleaf_set_t * restrict set);

int boxleaves_build_from_sorted_keys(const pkey_t * restrict sorted_keys,
                                     count_t N,
                                     count_t target,
                                     double tolerance,
                                     uint32_t max_side_cells,
                                     boxleaf_set_t * restrict out);

count_t boxleaves_find_leaf(const boxleaf_set_t * restrict set, pkey_t key);

void boxleaves_update_ranges_from_counts(boxleaf_set_t * restrict set,
                                         const count_t * restrict leaf_counts,
                                         const count_t * restrict gas_counts,
                                         const count_t * restrict offsets);

int boxleaves_verify_sorted_coverage(const boxleaf_set_t * restrict set,
                                     const pkey_t * restrict sorted_keys,
                                     count_t N);

#endif
