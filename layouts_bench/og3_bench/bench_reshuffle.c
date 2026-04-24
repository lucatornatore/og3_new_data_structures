/* =========================================================================
 * bench_reshuffle.c
 *
 * Global PH-order reshuffle benchmark. Dense mode builds a single global
 * permutation by sorting all particles by Peano-Hilbert key, with no type
 * grouping. Sparse mode builds a valid mostly-identity global permutation and
 * the layout implementations use changed-destination-only movement when the
 * permutation is sparse enough.
 * ========================================================================= */

#include "bench_common.h"
#include "layout_api.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void build_global_key_perm(const pkey_t * restrict keys,
                                  count_t N,
                                  count_t * restrict perm_out,
                                  pkey_t * restrict sorted_keys)
{
    count_t *iden = (count_t *)og3_aligned_alloc((size_t)N * sizeof(count_t));
    if (N && !iden) { fprintf(stderr, "build_global_key_perm: OOM\n"); exit(1); }
    for (count_t i = 0; i < N; ++i) iden[i] = i;

    /* radix_sort_u64_count is stable, so equal keys retain original-index
     * order via the iden[] input order. */
    radix_sort_u64_count(keys, iden, sorted_keys, perm_out, N);
    free(iden);
}

static int verify_global_sorted(const layout_ctx_t * restrict ctx, count_t N)
{
    for (count_t i = 1; i < N; ++i) {
        pkey_t a = layout_get_key(ctx, i - 1);
        pkey_t b = layout_get_key(ctx, i);
        if (a > b) {
            fprintf(stderr, "[verify] global key order failure at i=%llu (%llu > %llu)\n",
                    (unsigned long long)i,
                    (unsigned long long)a,
                    (unsigned long long)b);
            return 1;
        }
    }
    return 0;
}

static int verify_against_perm(const layout_ctx_t * restrict ctx,
                               const count_t * restrict perm,
                               const uint8_t * restrict types,
                               const ppid_t * restrict ids,
                               count_t N)
{
    for (count_t i = 0; i < N; ++i) {
        count_t old = perm[i];
        if (layout_get_id(ctx, i) != ids[old] || layout_get_type(ctx, i) != types[old]) {
            fprintf(stderr,
                    "[verify] semantic mismatch at slot %llu: got(id=%llu,type=%u) expected(id=%llu,type=%u)\n",
                    (unsigned long long)i,
                    (unsigned long long)layout_get_id(ctx, i),
                    (unsigned)layout_get_type(ctx, i),
                    (unsigned long long)ids[old],
                    (unsigned)types[old]);
            return 1;
        }
    }
    return 0;
}

static int exercise_key_pos_api(layout_ctx_t * restrict ctx, count_t N)
{
    if (N == 0) return 0;

    count_t probes[3] = {0, N / 2, N - 1};
    for (int p = 0; p < 3; ++p) {
        count_t i = probes[p];
        int already_done = 0;
        for (int q = 0; q < p; ++q) {
            if (probes[q] == i) already_done = 1;
        }
        if (already_done) continue;

        pkey_t old_key = layout_get_key(ctx, i);
        pkey_t new_key = old_key ^ UINT64_C(0x9e3779b97f4a7c15);
        if (new_key == old_key) new_key = old_key + 1u;

        layout_set_key(ctx, i, new_key);
        if (layout_get_key(ctx, i) != new_key) {
            fprintf(stderr, "[verify] layout_set_key failed at slot %llu\n",
                    (unsigned long long)i);
            return 1;
        }
        layout_set_key(ctx, i, old_key);
        if (layout_get_key(ctx, i) != old_key) {
            fprintf(stderr, "[verify] layout_set_key restore failed at slot %llu\n",
                    (unsigned long long)i);
            return 1;
        }

        pos_t old_pos[3];
        pos_t new_pos[3];
        pos_t got_pos[3];
        layout_get_pos(ctx, i, old_pos);
        new_pos[0] = old_pos[0] + (pos_t)1.25;
        new_pos[1] = old_pos[1] - (pos_t)2.50;
        new_pos[2] = old_pos[2] + (pos_t)3.75;

        layout_set_pos(ctx, i, new_pos);
        layout_get_pos(ctx, i, got_pos);
        if (got_pos[0] != new_pos[0] || got_pos[1] != new_pos[1] || got_pos[2] != new_pos[2]) {
            fprintf(stderr, "[verify] layout_set_pos failed at slot %llu\n",
                    (unsigned long long)i);
            return 1;
        }
        layout_set_pos(ctx, i, old_pos);
        layout_get_pos(ctx, i, got_pos);
        if (got_pos[0] != old_pos[0] || got_pos[1] != old_pos[1] || got_pos[2] != old_pos[2]) {
            fprintf(stderr, "[verify] layout_set_pos restore failed at slot %llu\n",
                    (unsigned long long)i);
            return 1;
        }
    }
    return 0;
}

int main(int argc, char **argv)
{
    count_t N        = (argc > 1) ? strtoull(argv[1], NULL, 10) : 1000000ULL;
    int     n_iter   = (argc > 2) ? atoi(argv[2])               : 10;
    double  fraction = (argc > 3) ? atof(argv[3])               : 1.0;
    uint32_t n_cells = 1u << 21;

    if (fraction <= 0.0 || fraction > 1.0) {
        fprintf(stderr, "fraction must be in (0,1]; got %g\n", fraction);
        return 1;
    }
    int dense = (fraction >= 0.9999);

    printf("\n================ bench_reshuffle  layout=%s ================\n", layout_name());
    printf("description : %s\n", layout_description());
    printf("N           : %llu\n", (unsigned long long)N);
    printf("n_iter      : %d\n", n_iter);
    printf("mode        : %s (fraction=%.4f)\n",
           dense ? "DENSE global PH sort" : "SPARSE global permutation", fraction);
    printf("usage       : bench_reshuffle_X [N] [n_iter] [fraction]\n\n");

    pcount_t c;
    counts_from_total(N, &c);
    printf("counts      : gas=%llu  dm=%llu  star=%llu  bh=%llu\n",
           (unsigned long long)c.n_gas, (unsigned long long)c.n_dm,
           (unsigned long long)c.n_star, (unsigned long long)c.n_bh);

    pkey_t *keys = (pkey_t *)og3_aligned_alloc((size_t)N * sizeof(pkey_t));
    pkey_t *sorted_keys = (pkey_t *)og3_aligned_alloc((size_t)N * sizeof(pkey_t));
    uint8_t *types = (uint8_t *)og3_aligned_alloc((size_t)N * sizeof(uint8_t));
    ppid_t *ids = (ppid_t *)og3_aligned_alloc((size_t)N * sizeof(ppid_t));
    count_t *perm = (count_t *)og3_aligned_alloc((size_t)N * sizeof(count_t));
    if (N && (!keys || !sorted_keys || !types || !ids || !perm)) {
        fprintf(stderr, "OOM\n");
        return 1;
    }

    gen_keys(keys, N, n_cells, 12345ULL);
    assign_types(types, &c, 67890ULL);
    for (count_t p = 0; p < N; ++p) ids[p] = (ppid_t)p;

    if (dense) {
        TIME_BLOCK("build perm (global PH sort)",
                   build_global_key_perm(keys, N, perm, sorted_keys));
        if (verify_sort(sorted_keys, N)) return 2;
    } else {
        TIME_BLOCK("build perm (sparse global)",
                   gen_sparse_permutation(perm, N, fraction, 98765ULL));
    }

    count_t changed = count_permutation_changed(perm, N);
    printf("perm changed slots      : %llu / %llu (%.4f)\n",
           (unsigned long long)changed,
           (unsigned long long)N,
           N ? (double)changed / (double)N : 0.0);

    layout_ctx_t *ctx = layout_alloc(&c);
    if (!ctx) { fprintf(stderr, "layout_alloc failed\n"); return 1; }
    TIME_BLOCK("layout_fill", layout_fill(ctx, keys, types, ids));
    if (exercise_key_pos_api(ctx, N)) {
        fprintf(stderr, "[verify] key/pos API self-test failed\n");
        layout_free(ctx);
        return 3;
    }

    TIME_BLOCK("warmup reshuffle", layout_reshuffle_full(ctx, perm));

    int failed = 0;
    failed |= verify_against_perm(ctx, perm, types, ids, N);
    if (dense) failed |= verify_global_sorted(ctx, N);
    failed |= layout_verify_deep(ctx);
    if (failed) {
        fprintf(stderr, "[verify] FAIL\n");
        layout_free(ctx);
        return 3;
    }
    printf("[verify] semantic alignment%s: OK\n", dense ? ", global PH order" : "");

    double t0 = cpu_time_p();
    for (int it = 0; it < n_iter; ++it)
        layout_reshuffle_full(ctx, perm);
    double t1 = cpu_time_p();

    double dt_total = t1 - t0;
    double ms_per_it = 1000.0 * dt_total / (double)n_iter;
    size_t dense_payload_bytes = layout_reshuffle_bytes(ctx);
    int streams = layout_reshuffle_streams(ctx);
    size_t index_upper = (size_t)streams * (size_t)N * sizeof(count_t);
    double mb_per_s = (ms_per_it > 0.0) ? (double)(dense_payload_bytes + index_upper) / ms_per_it * 1e-3 : 0.0;
    double ns_per_part = (N > 0) ? (double)ms_per_it * 1e6 / (double)N : 0.0;

    printf("\n---- results (%s, %s) ----\n", layout_name(), dense ? "dense" : "sparse");
    printf("streams permuted        : %d\n", streams);
    printf("payload copy bytes      : %zu  (%.1f MB, dense 4-pass estimate)\n",
           dense_payload_bytes, dense_payload_bytes * 1e-6);
    printf("perm-index bytes        : %zu  (%.1f MB, coarse upper estimate)\n",
           index_upper, index_upper * 1e-6);
    if (!dense) {
        double sparse_scale = N ? (double)changed / (2.0 * (double)N) : 0.0;
        printf("sparse payload estimate : %.1f MB (changed-destination copy path)\n",
               (double)dense_payload_bytes * sparse_scale * 1e-6);
    }
    printf("ms / reshuffle          : %.3f\n", ms_per_it);
    printf("ns / particle           : %.2f\n", ns_per_part);
    printf("throughput estimate     : %.2f MB/s  (%.2f GB/s)\n",
           mb_per_s, mb_per_s * 1e-3);

    layout_free(ctx);
    free(perm); free(keys); free(sorted_keys); free(types); free(ids);
    return 0;
}
