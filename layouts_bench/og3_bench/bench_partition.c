/* =========================================================================
 * bench_partition.c
 *
 * Same-type swap primitive benchmark. This remains a lower-level proxy for
 * future leaf-partition maintenance: it measures the cost of swapping
 * same-type particles without changing the type sequence at global slots.
 * ========================================================================= */

#include "bench_common.h"
#include "layout_api.h"

#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv)
{
    count_t N       = (argc > 1) ? strtoull(argv[1], NULL, 10) : 1000000ULL;
    double fraction = (argc > 2) ? atof(argv[2])               : 0.01;
    int n_iter      = (argc > 3) ? atoi(argv[3])               : 20;

    if (fraction <= 0.0 || fraction > 1.0) {
        fprintf(stderr, "fraction must be in (0, 1]\n");
        return 1;
    }

    printf("\n================ bench_partition  layout=%s ================\n", layout_name());
    printf("description : %s\n", layout_description());
    printf("N           : %llu\n", (unsigned long long)N);
    printf("fraction    : %.4f  (%llu swap pairs)\n", fraction,
           (unsigned long long)(fraction * (double)N));
    printf("n_iter      : %d\n\n", n_iter);

    pcount_t c;
    counts_from_total(N, &c);

    pkey_t *keys = (pkey_t *)og3_aligned_alloc((size_t)N * sizeof(pkey_t));
    uint8_t *types = (uint8_t *)og3_aligned_alloc((size_t)N * sizeof(uint8_t));
    ppid_t *ids = (ppid_t *)og3_aligned_alloc((size_t)N * sizeof(ppid_t));
    if (N && (!keys || !types || !ids)) {
        fprintf(stderr, "OOM input arrays\n");
        return 1;
    }

    gen_keys(keys, N, 1u << 21, 12345ULL);
    assign_types(types, &c, 67890ULL);
    for (count_t p = 0; p < N; ++p) ids[p] = (ppid_t)p;

    layout_ctx_t *ctx = layout_alloc(&c);
    if (!ctx) { fprintf(stderr, "layout_alloc failed\n"); return 1; }
    layout_fill(ctx, keys, types, ids);

    count_t n_pairs = (count_t)(fraction * (double)N);
    if (n_pairs == 0) { fprintf(stderr, "n_pairs == 0\n"); return 1; }
    swap_pair_t *pairs = (swap_pair_t *)og3_aligned_alloc((size_t)n_pairs * sizeof(swap_pair_t));
    if (!pairs) { fprintf(stderr, "OOM pairs\n"); return 1; }
    gen_swap_pairs_same_type(pairs, n_pairs, types, N, 24680ULL);

    size_t total_bytes_per_iter = 0;
    for (count_t i = 0; i < n_pairs; ++i) {
        uint8_t t = types[pairs[i].i];
        total_bytes_per_iter += layout_swap_bytes(ctx, t);
    }

    for (count_t i = 0; i < n_pairs; ++i)
        layout_swap_same_type(ctx, pairs[i].i, pairs[i].j);
    if (layout_verify_deep(ctx)) {
        fprintf(stderr, "[verify] FAIL after warmup swaps\n");
        return 2;
    }

    double t0 = cpu_time_p();
    for (int it = 0; it < n_iter; ++it) {
        for (count_t i = 0; i < n_pairs; ++i)
            layout_swap_same_type(ctx, pairs[i].i, pairs[i].j);
    }
    double t1 = cpu_time_p();

    double dt_total = t1 - t0;
    double ms_per_iter = 1000.0 * dt_total / (double)n_iter;
    double ns_per_swap = 1e9 * dt_total / ((double)n_iter * (double)n_pairs);
    double mb_per_s = (ms_per_iter > 0.0) ? (double)total_bytes_per_iter / ms_per_iter * 1e-3 : 0.0;

    printf("---- per-type swap cost ----\n");
    printf("  gas  : %zu B, %d streams\n", layout_swap_bytes(ctx, PT_GAS), layout_swap_streams(ctx, PT_GAS));
    printf("  dm   : %zu B, %d streams\n", layout_swap_bytes(ctx, PT_DM1), layout_swap_streams(ctx, PT_DM1));
    printf("  star : %zu B, %d streams\n", layout_swap_bytes(ctx, PT_STAR), layout_swap_streams(ctx, PT_STAR));
    printf("  bh   : %zu B, %d streams\n", layout_swap_bytes(ctx, PT_BH), layout_swap_streams(ctx, PT_BH));

    printf("\n---- results (%s) ----\n", layout_name());
    printf("pairs / iter            : %llu\n", (unsigned long long)n_pairs);
    printf("bytes / iter            : %zu  (%.1f MB)\n",
           total_bytes_per_iter, total_bytes_per_iter * 1e-6);
    printf("ms / iter               : %.3f\n", ms_per_iter);
    printf("ns / swap               : %.1f\n", ns_per_swap);
    printf("throughput              : %.2f MB/s  (%.2f GB/s)\n", mb_per_s, mb_per_s * 1e-3);

    free(pairs);
    layout_free(ctx);
    free(keys); free(types); free(ids);
    return 0;
}
