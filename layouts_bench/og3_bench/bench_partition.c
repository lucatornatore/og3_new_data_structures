#include "bench_common.h"
#include "layout_api.h"
#include "particle_init.h"
#include "boxleaves.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct { coord3u_t pos; pkey_t key; uint8_t type; ppid_t id; } seed_t;

static void build_perm_by_key(const pkey_t *keys, count_t N,
                              count_t *perm, pkey_t *sorted_keys)
{
    count_t *iden = (count_t *)og3_aligned_alloc((size_t)N * sizeof(count_t));
    if (N && !iden) { fprintf(stderr, "OOM iden\n"); exit(1); }
    for (count_t i = 0; i < N; ++i) iden[i] = i;
    radix_sort_u64_count(keys, iden, sorted_keys, perm, N);
    free(iden);
}

static uint32_t clamp_coord(int64_t v)
{
    if (v < 0) return 0;
    if ((uint64_t)v >= (uint64_t)OG3_PH_CELLS) return OG3_PH_CELLS - 1u;
    return (uint32_t)v;
}

static void perturb_one(layout_ctx_t *ctx, count_t i, uint32_t delta, rng_t *rng)
{
    pos_t p[3];
    layout_get_pos(ctx, i, p);
    int64_t span = (int64_t)delta;
    int64_t dx = span ? (int64_t)(rng_next(rng) % (uint64_t)(2 * span + 1)) - span : 0;
    int64_t dy = span ? (int64_t)(rng_next(rng) % (uint64_t)(2 * span + 1)) - span : 0;
    int64_t dz = span ? (int64_t)(rng_next(rng) % (uint64_t)(2 * span + 1)) - span : 0;
    uint32_t x = clamp_coord((int64_t)llround(p[0]) + dx);
    uint32_t y = clamp_coord((int64_t)llround(p[1]) + dy);
    uint32_t z = clamp_coord((int64_t)llround(p[2]) + dz);
    pos_t np[3] = {(pos_t)x, (pos_t)y, (pos_t)z};
    layout_set_pos(ctx, i, np);
    layout_set_key(ctx, i, ph_key_from_ijk(x, y, z));
}

static count_t build_partition_perm(const layout_ctx_t *ctx,
                                    boxleaf_set_t *leaves,
                                    count_t N,
                                    count_t *perm,
                                    count_t *leaf_of,
                                    count_t *leaf_counts,
                                    count_t *gas_counts,
                                    count_t *offsets,
                                    count_t *outside)
{
    const count_t nb = leaves->n_leaf;
    memset(leaf_counts, 0, (size_t)(nb + 1) * sizeof(count_t));
    memset(gas_counts, 0, (size_t)(nb + 1) * sizeof(count_t));
    *outside = 0;

    for (count_t i = 0; i < N; ++i) {
        count_t l = boxleaves_find_leaf(leaves, layout_get_key(ctx, i));
        if (l == BOXLEAF_INVALID) { l = nb; (*outside)++; }
        leaf_of[i] = l;
        leaf_counts[l]++;
        if (layout_get_type(ctx, i) == PT_GAS) gas_counts[l]++;
    }

    count_t run = 0;
    for (count_t l = 0; l <= nb; ++l) {
        offsets[l] = run;
        run += leaf_counts[l];
    }

    count_t *gas_cursor = (count_t *)og3_aligned_alloc((size_t)(nb + 1) * sizeof(count_t));
    count_t *other_cursor = (count_t *)og3_aligned_alloc((size_t)(nb + 1) * sizeof(count_t));
    if (!gas_cursor || !other_cursor) { fprintf(stderr, "OOM partition cursors\n"); exit(1); }

    for (count_t l = 0; l <= nb; ++l) {
        gas_cursor[l] = offsets[l];
        other_cursor[l] = offsets[l] + gas_counts[l];
    }

    /* Partition by leaf first and by gas/non-gas second. This keeps each leaf
     * contiguous and gives every leaf a gas subrange [gas_begin, gas_end). */
    for (count_t old = 0; old < N; ++old) {
        count_t l = leaf_of[old];
        if (layout_get_type(ctx, old) == PT_GAS)
            perm[gas_cursor[l]++] = old;
        else
            perm[other_cursor[l]++] = old;
    }

    boxleaves_update_ranges_from_counts(leaves, leaf_counts, gas_counts, offsets);

    free(gas_cursor);
    free(other_cursor);
    return count_permutation_changed(perm, N);
}

static int verify_partitioned(const layout_ctx_t *ctx,
                              const boxleaf_set_t *leaves,
                              count_t N,
                              count_t *outside)
{
    const count_t nb = leaves->n_leaf;
    count_t prev = 0;
    *outside = 0;

    for (count_t i = 0; i < N; ++i) {
        count_t l = boxleaves_find_leaf(leaves, layout_get_key(ctx, i));
        if (l == BOXLEAF_INVALID) { l = nb; (*outside)++; }
        if (i && l < prev) {
            fprintf(stderr, "[verify] non-contiguous leaf order at %llu\n",
                    (unsigned long long)i);
            return 1;
        }
        prev = l;
    }

    for (count_t l = 0; l < nb; ++l) {
        const boxleaf_t *b = &leaves->leaves[l];
        if (b->begin > b->gas_begin || b->gas_begin > b->gas_end || b->gas_end > b->end) {
            fprintf(stderr, "[verify] invalid gas range in leaf %llu\n",
                    (unsigned long long)l);
            return 1;
        }
        for (count_t i = b->begin; i < b->end; ++i) {
            count_t li = boxleaves_find_leaf(leaves, layout_get_key(ctx, i));
            if (li != l) {
                fprintf(stderr, "[verify] particle %llu in leaf range %llu maps to %llu\n",
                        (unsigned long long)i, (unsigned long long)l, (unsigned long long)li);
                return 1;
            }
            uint8_t t = layout_get_type(ctx, i);
            if (i < b->gas_end && t != PT_GAS) {
                fprintf(stderr, "[verify] non-gas in gas subrange at %llu\n",
                        (unsigned long long)i);
                return 1;
            }
            if (i >= b->gas_end && t == PT_GAS) {
                fprintf(stderr, "[verify] gas after gas subrange at %llu\n",
                        (unsigned long long)i);
                return 1;
            }
        }
    }

    return 0;
}

int main(int argc, char **argv)
{
    count_t N = (argc > 1) ? strtoull(argv[1], NULL, 10) : 1000000ULL;
    int n_iter = (argc > 2) ? atoi(argv[2]) : 5;
    uint32_t perturb = (argc > 3) ? (uint32_t)strtoul(argv[3], NULL, 10) : 256u;
    count_t target = (argc > 4) ? strtoull(argv[4], NULL, 10) : 32ULL;
    double tol = (argc > 5) ? atof(argv[5]) : 0.5;
    uint32_t max_side = (argc > 6) ? (uint32_t)strtoul(argv[6], NULL, 10) : (1u << OG3_PH_BITS);
    int mode_i = (argc > 7) ? atoi(argv[7]) : 0;
    double temp = (argc > 8) ? atof(argv[8]) : 0.25;
    og3_pos_mode_t mode = mode_i ? OG3_POS_CLUSTERED : OG3_POS_PLAIN;

    printf("\n================ bench_partition  layout=%s ================\n", layout_name());
    printf("description : %s\n", layout_description());
    printf("N=%llu n_iter=%d perturb=+/- %u target=%llu tol=%.3f max_side=%u mode=%s temperature=%.4f\n",
           (unsigned long long)N, n_iter, perturb, (unsigned long long)target, tol, max_side,
           mode == OG3_POS_CLUSTERED ? "clustered" : "plain", temp);
    printf("usage: bench_partition_X [N] [n_iter] [perturb_cells] [target_leaf] [tol] [max_side] [mode:0|1] [temperature]\n\n");

    pcount_t c;
    counts_from_total(N, &c);

    coord3u_t *pos = (coord3u_t *)og3_aligned_alloc((size_t)N * sizeof(coord3u_t));
    pkey_t *keys = (pkey_t *)og3_aligned_alloc((size_t)N * sizeof(pkey_t));
    pkey_t *sorted_keys = (pkey_t *)og3_aligned_alloc((size_t)N * sizeof(pkey_t));
    uint8_t *types = (uint8_t *)og3_aligned_alloc((size_t)N * sizeof(uint8_t));
    ppid_t *ids = (ppid_t *)og3_aligned_alloc((size_t)N * sizeof(ppid_t));
    count_t *sort_perm = (count_t *)og3_aligned_alloc((size_t)N * sizeof(count_t));
    seed_t *sorted = (seed_t *)og3_aligned_alloc((size_t)N * sizeof(seed_t));
    if (N && (!pos || !keys || !sorted_keys || !types || !ids || !sort_perm || !sorted)) {
        fprintf(stderr, "OOM inputs\n");
        return 1;
    }

    TIME_BLOCK("generate positions", og3_generate_positions(pos, N, mode, temp, 12345ULL));
    TIME_BLOCK("derive PH keys", og3_keys_from_positions(pos, keys, N));
    assign_types(types, &c, 67890ULL);
    for (count_t i = 0; i < N; ++i) ids[i] = (ppid_t)i;

    TIME_BLOCK("global PH sort", build_perm_by_key(keys, N, sort_perm, sorted_keys));
    if (verify_sort(sorted_keys, N)) return 2;
    for (count_t i = 0; i < N; ++i) {
        count_t o = sort_perm[i];
        sorted[i] = (seed_t){pos[o], keys[o], types[o], ids[o]};
    }

    boxleaf_set_t leaves;
    TIME_BLOCK("build box leaves", boxleaves_build_from_sorted_keys(sorted_keys, N, target, tol, max_side, &leaves));
    printf("box leaves: %llu target=%llu min=%llu max=%llu underfull=%llu overfull=%llu\n",
           (unsigned long long)leaves.n_leaf,
           (unsigned long long)leaves.target,
           (unsigned long long)leaves.min_particles,
           (unsigned long long)leaves.max_particles,
           (unsigned long long)leaves.underfull,
           (unsigned long long)leaves.overfull);

    for (count_t i = 0; i < N; ++i) {
        pos[i] = sorted[i].pos;
        keys[i] = sorted[i].key;
        types[i] = sorted[i].type;
        ids[i] = sorted[i].id;
    }

    layout_ctx_t *ctx = layout_alloc(&c);
    if (!ctx) { fprintf(stderr, "layout_alloc failed\n"); return 1; }
    TIME_BLOCK("layout_fill sorted", layout_fill(ctx, keys, types, ids));
    og3_set_layout_positions(ctx, pos, N);
    if (layout_verify_deep(ctx)) { fprintf(stderr, "initial deep verify failed\n"); return 3; }

    count_t *perm = (count_t *)og3_aligned_alloc((size_t)N * sizeof(count_t));
    count_t *leaf_of = (count_t *)og3_aligned_alloc((size_t)N * sizeof(count_t));
    count_t *leaf_counts = (count_t *)og3_aligned_alloc((size_t)(leaves.n_leaf + 1) * sizeof(count_t));
    count_t *gas_counts = (count_t *)og3_aligned_alloc((size_t)(leaves.n_leaf + 1) * sizeof(count_t));
    count_t *offsets = (count_t *)og3_aligned_alloc((size_t)(leaves.n_leaf + 1) * sizeof(count_t));
    if (N && (!perm || !leaf_of || !leaf_counts || !gas_counts || !offsets)) {
        fprintf(stderr, "OOM partition\n");
        return 1;
    }

    rng_t rng;
    rng_seed(&rng, 24680ULL);
    count_t changed = 0, outside = 0;
    double t0 = cpu_time_p();
    for (int it = 0; it < n_iter; ++it) {
        for (count_t i = 0; i < N; ++i) perturb_one(ctx, i, perturb, &rng);
        changed = build_partition_perm(ctx, &leaves, N, perm, leaf_of,
                                       leaf_counts, gas_counts, offsets, &outside);
        layout_reshuffle_full(ctx, perm);
    }
    double t1 = cpu_time_p();

    count_t outside_v = 0;
    if (verify_partitioned(ctx, &leaves, N, &outside_v) || layout_verify_deep(ctx)) {
        fprintf(stderr, "[verify] FAIL\n");
        return 4;
    }
    printf("[verify] leaf partitioning + gas subranges: OK\n");

    double ms = n_iter ? 1000.0 * (t1 - t0) / (double)n_iter : 0.0;
    size_t bytes = layout_reshuffle_bytes(ctx);
    size_t idx_bytes = (size_t)layout_reshuffle_streams(ctx) * (size_t)N * sizeof(count_t);
    printf("\n---- results (%s) ----\n", layout_name());
    printf("last changed slots: %llu / %llu (%.4f)\n",
           (unsigned long long)changed, (unsigned long long)N,
           N ? (double)changed / (double)N : 0.0);
    printf("outside old leaves: %llu last, %llu verify\n",
           (unsigned long long)outside, (unsigned long long)outside_v);
    printf("payload copy bytes: %zu (%.1f MB)\n", bytes, bytes * 1e-6);
    printf("perm-index bytes  : %zu (%.1f MB)\n", idx_bytes, idx_bytes * 1e-6);
    printf("ms/iter           : %.3f\n", ms);
    printf("ns/particle       : %.2f\n", N ? ms * 1e6 / (double)N : 0.0);

    free(perm);
    free(leaf_of);
    free(leaf_counts);
    free(gas_counts);
    free(offsets);
    layout_free(ctx);
    boxleaf_set_free(&leaves);
    free(sorted);
    free(sort_perm);
    free(pos);
    free(keys);
    free(sorted_keys);
    free(types);
    free(ids);
    return 0;
}
