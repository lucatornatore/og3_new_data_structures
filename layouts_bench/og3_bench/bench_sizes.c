/* =========================================================================
 * bench_sizes.c
 *
 * Prints actual sizeof/alignment-dependent accounting for the selected layout.
 * Build creates bench_sizes_{A,Bp,Bc,C}, each linked to one layout object.
 * ========================================================================= */

#include "bench_common.h"
#include "fields.h"
#include "layout_api.h"

#include <stdio.h>
#include <stdlib.h>

static void print_shared_atoms(void)
{
    printf("OG3_ATOM_ALIGN = %d\n", OG3_ATOM_ALIGN);
    printf("\nshared fine atoms:\n");
    printf("  PCore    %5zu\n", sizeof(PCore));
    printf("  PDyn     %5zu\n", sizeof(PDyn));
    printf("  PTime    %5zu\n", sizeof(PTime));
    printf("  PMeta    %5zu\n", sizeof(PMeta));
    printf("  PLinkage %5zu\n", sizeof(PLinkage));
    printf("  GasCore  %5zu\n", sizeof(GasCore));
    printf("  GasGrad  %5zu\n", sizeof(GasGrad));
    printf("  GasMetal %5zu\n", sizeof(GasMetal));
    printf("  GasSF    %5zu\n", sizeof(GasSF));
    printf("  StarCore %5zu\n", sizeof(StarCore));
    printf("  StarMeta %5zu\n", sizeof(StarMeta));
    printf("  BHCore   %5zu\n", sizeof(BHCore));
    printf("  BHEnv    %5zu\n", sizeof(BHEnv));
    printf("  BHRepos  %5zu\n", sizeof(BHRepos));
}

int main(int argc, char **argv)
{
    count_t N = (argc > 1) ? strtoull(argv[1], NULL, 10) : 1000000ULL;
    pcount_t c;
    counts_from_total(N, &c);

    printf("\n================ bench_sizes layout=%s ================\n", layout_name());
    printf("description : %s\n", layout_description());
    printf("N           : %llu\n", (unsigned long long)N);
    printf("counts      : gas=%llu  dm=%llu  star=%llu  bh=%llu\n",
           (unsigned long long)c.n_gas, (unsigned long long)c.n_dm,
           (unsigned long long)c.n_star, (unsigned long long)c.n_bh);
    print_shared_atoms();

    layout_ctx_t *ctx = layout_alloc(&c);
    if (!ctx) {
        fprintf(stderr, "layout_alloc failed\n");
        return 1;
    }

    printf("\nlayout accounting:\n");
    printf("  reshuffle streams          %d\n", layout_reshuffle_streams(ctx));
    printf("  dense payload-copy bytes   %zu (%.3f MB)\n",
           layout_reshuffle_bytes(ctx), layout_reshuffle_bytes(ctx) * 1e-6);
    printf("  swap gas bytes             %zu\n", layout_swap_bytes(ctx, PT_GAS));
    printf("  swap dm bytes              %zu\n", layout_swap_bytes(ctx, PT_DM1));
    printf("  swap star bytes            %zu\n", layout_swap_bytes(ctx, PT_STAR));
    printf("  swap bh bytes              %zu\n", layout_swap_bytes(ctx, PT_BH));

    layout_free(ctx);
    return 0;
}
