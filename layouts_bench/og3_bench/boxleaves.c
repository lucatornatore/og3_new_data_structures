#include "boxleaves.h"

#include <string.h>

/* -------------------------------------------------------------------------
 * Box-leaf builder.
 *
 * This module keeps the code standalone but follows the logic of the uploaded
 * create_boxleaves prototype:
 *   - keys are already globally PH-sorted;
 *   - leaves are PH cubes;
 *   - splitting is driven by target/min/max occupancy plus a maximum allowed
 *     spatial side length;
 *   - every emitted leaf owns a contiguous particle interval.
 *
 * Two deliberate corrections relative to the uploaded prototype:
 *   1. small-N allocation cannot underflow/allocate zero leaves;
 *   2. boundary cases never read K[pstart+1] unless pstart+1 < N.
 * ------------------------------------------------------------------------- */

static count_t ceil_count_from_double(double x)
{
    if (x <= 1.0) return 1;
    count_t y = (count_t)x;
    return ((double)y < x) ? y + 1 : y;
}

static count_t floor_count_from_double(double x)
{
    if (x <= 1.0) return 1;
    return (count_t)x;
}

static count_t min_count(count_t target, double tol)
{
    if (target == 0) target = 1;
    if (tol < 0.0) tol = 0.0;
    if (tol > 1.0) tol = 1.0;
    return floor_count_from_double((double)target * (1.0 - tol));
}

static count_t max_count(count_t target, double tol)
{
    if (target == 0) target = 1;
    if (tol < 0.0) tol = 0.0;
    return ceil_count_from_double((double)target * (1.0 + tol));
}

static void reserve_leaf(boxleaf_set_t *set, count_t extra)
{
    if (set->n_leaf + extra <= set->capacity) return;
    count_t nc = set->capacity ? set->capacity * 2 : 1024;
    while (nc < set->n_leaf + extra) nc *= 2;
    boxleaf_t *p = (boxleaf_t *)realloc(set->leaves, (size_t)nc * sizeof(boxleaf_t));
    if (!p) { fprintf(stderr, "boxleaves: OOM while reserving leaves\n"); exit(1); }
    set->leaves = p;
    set->capacity = nc;
}

void boxleaf_set_free(boxleaf_set_t * restrict set)
{
    if (!set) return;
    free(set->leaves);
    memset(set, 0, sizeof(*set));
}

static uint32_t side_for_level(unsigned level)
{
    if (level >= OG3_PH_BITS) return 1u;
    return 1u << (OG3_PH_BITS - level);
}

static pkey_t key_begin_from_prefix(uint64_t prefix, unsigned level)
{
    if (level == 0) return 0;
    return (pkey_t)(prefix << (OG3_PH_KEY_BITS - 3u * level));
}

static pkey_t key_end_from_prefix(uint64_t prefix, unsigned level)
{
    if (level == 0) return (UINT64_C(1) << OG3_PH_KEY_BITS);
    return (pkey_t)((prefix + 1u) << (OG3_PH_KEY_BITS - 3u * level));
}

static uint64_t prefix_from_key(pkey_t key, unsigned level)
{
    if (level == 0) return 0;
    return (uint64_t)(key >> (OG3_PH_KEY_BITS - 3u * level));
}

static unsigned max_common_level(pkey_t a, pkey_t b)
{
    pkey_t x = a ^ b;
    if (x == 0) return OG3_PH_BITS;
    unsigned used_lz = (unsigned)__builtin_clzll(x) - (64u - OG3_PH_KEY_BITS);
    return used_lz / 3u;
}

static unsigned level_for_side(uint32_t max_side_cells)
{
    if (max_side_cells == 0) max_side_cells = 1;
    if (max_side_cells >= (1u << OG3_PH_BITS)) return 0;
    unsigned level = 0;
    uint32_t side = 1u << OG3_PH_BITS;
    while (level < OG3_PH_BITS && side > max_side_cells) {
        side >>= 1;
        level++;
    }
    return level;
}

static void emit_leaf(boxleaf_set_t *set,
                      const pkey_t *keys,
                      count_t begin,
                      count_t end,
                      unsigned level)
{
    if (begin >= end) return;
    uint64_t prefix = prefix_from_key(keys[begin], level);
    pkey_t kb = key_begin_from_prefix(prefix, level);
    pkey_t ke = key_end_from_prefix(prefix, level);

    reserve_leaf(set, 1);
    boxleaf_t *l = &set->leaves[set->n_leaf++];
    l->key_begin = kb;
    l->key_end = ke;
    l->begin = begin;
    l->end = end;
    l->gas_begin = begin;
    l->gas_end = begin;
    l->octant = kb;
    l->side_cells = side_for_level(level);
    l->level = (uint8_t)level;

    count_t n = end - begin;
    if (n < set->min_particles) set->underfull++;
    if (n > set->max_particles) set->overfull++;
    if (l->side_cells <= set->max_side_cells && n < set->min_particles)
        set->spatially_forced += 0; /* kept for accounting symmetry */
}

static count_t first_not_in_prefix(const pkey_t *keys,
                                   count_t begin,
                                   count_t end,
                                   unsigned level)
{
    uint64_t prefix = prefix_from_key(keys[begin], level);
    count_t p = begin + 1;
    while (p < end && prefix_from_key(keys[p], level) == prefix) p++;
    return p;
}

static void build_range_user_style(const pkey_t *keys,
                                   count_t begin,
                                   count_t end,
                                   unsigned min_level,
                                   boxleaf_set_t *set)
{
    count_t pstart = begin;
    const unsigned side_level = level_for_side(set->max_side_cells);

    while (pstart < end) {
        count_t remaining = end - pstart;
        count_t guess_n = remaining < set->target ? remaining : set->target;
        if (guess_n == 0) guess_n = 1;

        count_t plast = pstart + guess_n - 1;
        if (plast >= end) plast = end - 1;

        unsigned level = max_common_level(keys[pstart], keys[plast]);
        if (level < min_level) level = min_level;
        if (level < side_level) level = side_level;
        if (level > OG3_PH_BITS) level = OG3_PH_BITS;

        count_t pend = first_not_in_prefix(keys, pstart, end, level);
        count_t bunch = pend - pstart;

        /* If this PH cube is overpopulated, refine it. The uploaded code used
         * min/max particle thresholds; here tolerance maps target -> [min,max].
         * We never split below one integer cell per axis. */
        while (bunch > set->max_particles && level < OG3_PH_BITS) {
            unsigned next_level = level + 1;
            count_t next_pend = first_not_in_prefix(keys, pstart, end, next_level);
            count_t next_bunch = next_pend - pstart;

            level = next_level;
            pend = next_pend;
            bunch = next_bunch;
        }

        emit_leaf(set, keys, pstart, pend, level);

        if (pend <= pstart) { /* should be impossible, but prevents dead loops */
            fprintf(stderr, "boxleaves: internal non-progress at particle %llu\n",
                    (unsigned long long)pstart);
            exit(2);
        }

        if (pend < end) {
            unsigned common = max_common_level(keys[pstart], keys[pend]);
            min_level = common + 1u;
            if (min_level > OG3_PH_BITS) min_level = OG3_PH_BITS;
            if (min_level < side_level) min_level = side_level;
        }
        pstart = pend;
    }
}

int boxleaves_build_from_sorted_keys(const pkey_t * restrict sorted_keys,
                                     count_t N,
                                     count_t target,
                                     double tolerance,
                                     uint32_t max_side_cells,
                                     boxleaf_set_t * restrict out)
{
    memset(out, 0, sizeof(*out));
    if (target == 0) target = 1;
    if (tolerance < 0.0) tolerance = 0.0;
    if (tolerance > 1.0) tolerance = 1.0;
    if (max_side_cells == 0) max_side_cells = 1u;
    if (max_side_cells > (1u << OG3_PH_BITS)) max_side_cells = (1u << OG3_PH_BITS);

    out->target = target;
    out->tolerance = tolerance;
    out->min_particles = min_count(target, tolerance);
    out->max_particles = max_count(target, tolerance);
    out->max_side_cells = max_side_cells;

    if (N == 0) return 0;

    /* Small-N fix: capacity must be at least one even when N << target. */
    out->capacity = (N / target) * 2u + 8u;
    if (out->capacity < 8u) out->capacity = 8u;
    out->leaves = (boxleaf_t *)malloc((size_t)out->capacity * sizeof(boxleaf_t));
    if (!out->leaves) { fprintf(stderr, "boxleaves: OOM initial allocation\n"); return 1; }

    build_range_user_style(sorted_keys, 0, N, 0, out);
    return boxleaves_verify_sorted_coverage(out, sorted_keys, N);
}

count_t boxleaves_find_leaf(const boxleaf_set_t * restrict set, pkey_t key)
{
    count_t lo = 0, hi = set->n_leaf;
    while (lo < hi) {
        count_t mid = lo + (hi - lo) / 2u;
        if (key < set->leaves[mid].key_begin) hi = mid;
        else if (key >= set->leaves[mid].key_end) lo = mid + 1u;
        else return mid;
    }
    return BOXLEAF_INVALID;
}

void boxleaves_update_ranges_from_counts(boxleaf_set_t * restrict set,
                                         const count_t * restrict leaf_counts,
                                         const count_t * restrict gas_counts,
                                         const count_t * restrict offsets)
{
    for (count_t l = 0; l < set->n_leaf; ++l) {
        boxleaf_t *b = &set->leaves[l];
        b->begin = offsets[l];
        b->end = offsets[l] + leaf_counts[l];
        b->gas_begin = offsets[l];
        b->gas_end = offsets[l] + gas_counts[l];
    }
}

int boxleaves_verify_sorted_coverage(const boxleaf_set_t * restrict set,
                                     const pkey_t * restrict sorted_keys,
                                     count_t N)
{
    if (N == 0) return 0;
    if (!set || set->n_leaf == 0) {
        fprintf(stderr, "boxleaves: no leaves for non-empty particle set\n");
        return 1;
    }

    count_t prev_end = 0;
    pkey_t prev_key_end = 0;
    for (count_t l = 0; l < set->n_leaf; ++l) {
        const boxleaf_t *b = &set->leaves[l];
        if (b->begin != prev_end || b->begin >= b->end || b->end > N) {
            fprintf(stderr, "boxleaves: particle range coverage failure at leaf %llu\n",
                    (unsigned long long)l);
            return 1;
        }
        if (b->key_begin >= b->key_end) {
            fprintf(stderr, "boxleaves: invalid key interval at leaf %llu\n",
                    (unsigned long long)l);
            return 1;
        }
        if (sorted_keys[b->begin] < b->key_begin || sorted_keys[b->end - 1] >= b->key_end) {
            fprintf(stderr, "boxleaves: key coverage failure at leaf %llu\n",
                    (unsigned long long)l);
            return 1;
        }
        if (l > 0 && b->key_begin < prev_key_end) {
            fprintf(stderr, "boxleaves: overlapping PH intervals at leaf %llu\n",
                    (unsigned long long)l);
            return 1;
        }
        if (b->side_cells > set->max_side_cells) {
            fprintf(stderr, "boxleaves: max-side failure at leaf %llu: side=%u max=%u\n",
                    (unsigned long long)l, b->side_cells, set->max_side_cells);
            return 1;
        }
        prev_end = b->end;
        prev_key_end = b->key_end;
    }
    if (prev_end != N) {
        fprintf(stderr, "boxleaves: final particle coverage failure\n");
        return 1;
    }
    return 0;
}
