/* =========================================================================
 * layout_Bc.c
 *
 * Layout Bc: coarse common atoms + one composite type-specific struct per
 * type + explicit cross-index. Global PH reshuffle permutes four common
 * streams; type-specific composites stay put.
 * ========================================================================= */

#include "layout_api.h"
#include "fields.h"
#include "bench_common.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct OG3_ATOM {
    pos_t   pos[3];
    real_t  mass;
    pkey_t  key;
    uint8_t type;
    uint8_t flags;
    int16_t time_bin;
    int32_t leaf_idx;
} PCoreB;

typedef struct OG3_ATOM {
    real_t  vel[3];
    real_t  dp[3];
    acc_t   grav_acc[3];
    acc_t   grav_pm[3];
    real_t  old_acc;
    real_t  num_ngb;
    int32_t true_ngb;
} PDynB;

typedef struct OG3_ATOM {
    ppid_t     id;
    time_int_t ti_begin;
    time_int_t ti_current;
    int32_t    dt_step;
    real_t     grav_cost[6];
} PAuxB;

typedef struct OG3_ATOM {
    count_t type_idx;
} PLinkageB;

typedef struct OG3_ATOM {
    GasCore  core;
    GasGrad  grad;
    GasMetal metal;
    GasSF    sf;
} GasAllB;

typedef struct OG3_ATOM {
    StarCore core;
    StarMeta meta;
} StarAllB;

typedef struct OG3_ATOM {
    BHCore  core;
    BHEnv   env;
    BHRepos repos;
} BHAllB;

struct layout_ctx {
    pcount_t counts;

    PCoreB    *core;
    PDynB     *dyn;
    PAuxB     *aux;
    PLinkageB *linkage;

    GasAllB  *gas_all;
    StarAllB *star_all;
    BHAllB   *bh_all;

    void     *scratch;
    size_t    scratch_bytes;
    count_t  *scratch_idx;
};

const char *layout_name(void)        { return "Bc"; }
const char *layout_description(void) { return "Coarse common + composite per-type + cross-index"; }

layout_ctx_t *layout_alloc(const pcount_t * restrict c)
{
    layout_ctx_t *ctx = (layout_ctx_t *)calloc(1, sizeof(*ctx));
    if (!ctx) return NULL;
    ctx->counts = *c;
    count_t N = total_parts(c);

    ctx->core    = (PCoreB    *)og3_aligned_alloc((size_t)N * sizeof(PCoreB));
    ctx->dyn     = (PDynB     *)og3_aligned_alloc((size_t)N * sizeof(PDynB));
    ctx->aux     = (PAuxB     *)og3_aligned_alloc((size_t)N * sizeof(PAuxB));
    ctx->linkage = (PLinkageB *)og3_aligned_alloc((size_t)N * sizeof(PLinkageB));

    if (c->n_gas)  ctx->gas_all  = (GasAllB  *)og3_aligned_alloc((size_t)c->n_gas * sizeof(GasAllB));
    if (c->n_star) ctx->star_all = (StarAllB *)og3_aligned_alloc((size_t)c->n_star * sizeof(StarAllB));
    if (c->n_bh)   ctx->bh_all   = (BHAllB   *)og3_aligned_alloc((size_t)c->n_bh * sizeof(BHAllB));

    /* Scratch only needs to hold one common stream. Type-specific composites
     * are not permuted by the global reshuffle in Bc. */
    size_t max_common = sizeof(PCoreB);
    if (sizeof(PDynB)     > max_common) max_common = sizeof(PDynB);
    if (sizeof(PAuxB)     > max_common) max_common = sizeof(PAuxB);
    if (sizeof(PLinkageB) > max_common) max_common = sizeof(PLinkageB);
    ctx->scratch_bytes = (size_t)N * max_common;
    ctx->scratch = og3_aligned_alloc(ctx->scratch_bytes);
    ctx->scratch_idx = (count_t *)og3_aligned_alloc((size_t)N * sizeof(count_t));

    if ((N && (!ctx->core || !ctx->dyn || !ctx->aux || !ctx->linkage ||
               !ctx->scratch || !ctx->scratch_idx)) ||
        (c->n_gas && !ctx->gas_all) ||
        (c->n_star && !ctx->star_all) ||
        (c->n_bh && !ctx->bh_all)) {
        layout_free(ctx);
        return NULL;
    }
    return ctx;
}

void layout_free(layout_ctx_t * restrict ctx)
{
    if (!ctx) return;
    free(ctx->core); free(ctx->dyn); free(ctx->aux); free(ctx->linkage);
    free(ctx->gas_all); free(ctx->star_all); free(ctx->bh_all);
    free(ctx->scratch); free(ctx->scratch_idx);
    free(ctx);
}

void layout_fill(layout_ctx_t * restrict ctx,
                 const pkey_t * restrict keys,
                 const uint8_t * restrict types,
                 const ppid_t * restrict ids)
{
    count_t N = total_parts(&ctx->counts);
    count_t ig = 0, is = 0, ib = 0;

    for (count_t p = 0; p < N; ++p) {
        uint8_t t = types[p];

        memset(&ctx->core[p], 0, sizeof(PCoreB));
        ctx->core[p].key = keys[p];
        ctx->core[p].type = t;
        ctx->core[p].pos[0] = (pos_t)(p * 1e-6);

        memset(&ctx->dyn[p], 0, sizeof(PDynB));
        ctx->dyn[p].vel[0] = (real_t)(p * 1e-3);

        memset(&ctx->aux[p], 0, sizeof(PAuxB));
        ctx->aux[p].id = ids[p];
        ctx->aux[p].ti_current = (time_int_t)p;

        memset(&ctx->linkage[p], 0, sizeof(PLinkageB));
        switch (t) {
        case PT_GAS:
            ctx->linkage[p].type_idx = ig;
            memset(&ctx->gas_all[ig], 0, sizeof(GasAllB));
            ctx->gas_all[ig].core.density = (real_t)(p * 1e-4);
            ctx->gas_all[ig].metal.metals[0] = (real_t)(p % 1000);
#ifdef DEBUG
            ctx->gas_all[ig].sf.swallow_id = og3_debug_marker64(ids[p], 1u);
#endif
            ++ig;
            break;
        case PT_STAR:
            ctx->linkage[p].type_idx = is;
            memset(&ctx->star_all[is], 0, sizeof(StarAllB));
            ctx->star_all[is].core.stellar_age = (real_t)(p * 1e-5);
#ifdef DEBUG
            ctx->star_all[is].core.pid = og3_debug_marker32(ids[p], 2u);
#endif
            ++is;
            break;
        case PT_BH:
            ctx->linkage[p].type_idx = ib;
            memset(&ctx->bh_all[ib], 0, sizeof(BHAllB));
            ctx->bh_all[ib].core.bh_mass = (double)(p * 1e-6);
#ifdef DEBUG
            ctx->bh_all[ib].core.swallow_id = og3_debug_marker64(ids[p], 3u);
#endif
            ++ib;
            break;
        default:
            ctx->linkage[p].type_idx = 0;
            break;
        }
    }
    assert(ig == ctx->counts.n_gas);
    assert(is == ctx->counts.n_star);
    assert(ib == ctx->counts.n_bh);
}

void layout_reshuffle_full(layout_ctx_t * restrict ctx, const count_t * restrict perm)
{
    count_t N = total_parts(&ctx->counts);
    og3_apply_perm_bytes(ctx->core,    perm, N, sizeof(PCoreB),    ctx->scratch, ctx->scratch_idx);
    og3_apply_perm_bytes(ctx->dyn,     perm, N, sizeof(PDynB),     ctx->scratch, ctx->scratch_idx);
    og3_apply_perm_bytes(ctx->aux,     perm, N, sizeof(PAuxB),     ctx->scratch, ctx->scratch_idx);
    og3_apply_perm_bytes(ctx->linkage, perm, N, sizeof(PLinkageB), ctx->scratch, ctx->scratch_idx);
}

#define SWAP_STRUCT(a, b, S) do {    \
    S tmp;                           \
    memcpy(&tmp, &(a), sizeof(S));   \
    memcpy(&(a), &(b), sizeof(S));   \
    memcpy(&(b), &tmp, sizeof(S));   \
} while (0)

void layout_swap_same_type(layout_ctx_t * restrict ctx, count_t i, count_t j)
{
    SWAP_STRUCT(ctx->core[i],    ctx->core[j],    PCoreB);
    SWAP_STRUCT(ctx->dyn[i],     ctx->dyn[j],     PDynB);
    SWAP_STRUCT(ctx->aux[i],     ctx->aux[j],     PAuxB);
    SWAP_STRUCT(ctx->linkage[i], ctx->linkage[j], PLinkageB);
}

size_t layout_reshuffle_bytes(const layout_ctx_t * restrict ctx)
{
    count_t N = total_parts(&ctx->counts);
    size_t per = sizeof(PCoreB) + sizeof(PDynB) + sizeof(PAuxB) + sizeof(PLinkageB);
    return og3_perm_dense_copy_bytes(N, per);
}

int layout_reshuffle_streams(const layout_ctx_t * restrict ctx)
{
    (void)ctx;
    return 4;
}

size_t layout_swap_bytes(const layout_ctx_t * restrict ctx, uint8_t t)
{
    (void)ctx; (void)t;
    return 6ull * (sizeof(PCoreB) + sizeof(PDynB) + sizeof(PAuxB) + sizeof(PLinkageB));
}

int layout_swap_streams(const layout_ctx_t * restrict ctx, uint8_t t)
{
    (void)ctx; (void)t;
    return 4;
}

pkey_t layout_get_key(const layout_ctx_t * restrict ctx, count_t i) { return ctx->core[i].key; }

void layout_set_key(layout_ctx_t * restrict ctx, count_t i, pkey_t key)
{
    ctx->core[i].key = key;
}

void layout_get_pos(const layout_ctx_t * restrict ctx, count_t i,
                    pos_t out[restrict 3])
{
    out[0] = ctx->core[i].pos[0];
    out[1] = ctx->core[i].pos[1];
    out[2] = ctx->core[i].pos[2];
}

void layout_set_pos(layout_ctx_t * restrict ctx, count_t i,
                    const pos_t in[restrict 3])
{
    ctx->core[i].pos[0] = in[0];
    ctx->core[i].pos[1] = in[1];
    ctx->core[i].pos[2] = in[2];
}

uint8_t layout_get_type(const layout_ctx_t * restrict ctx, count_t i) { return ctx->core[i].type; }
ppid_t layout_get_id(const layout_ctx_t * restrict ctx, count_t i) { return ctx->aux[i].id; }

int layout_verify_deep(const layout_ctx_t * restrict ctx)
{
#ifndef DEBUG
    (void)ctx;
    return 0;
#else
    count_t N = total_parts(&ctx->counts);
    int errors = 0;
    for (count_t i = 0; i < N; ++i) {
        uint8_t t = ctx->core[i].type;
        ppid_t id = ctx->aux[i].id;
        count_t k = ctx->linkage[i].type_idx;
        if (t == PT_GAS) {
            if (k >= ctx->counts.n_gas || ctx->gas_all[k].sf.swallow_id != og3_debug_marker64(id, 1u)) errors++;
        } else if (t == PT_STAR) {
            if (k >= ctx->counts.n_star || ctx->star_all[k].core.pid != og3_debug_marker32(id, 2u)) errors++;
        } else if (t == PT_BH) {
            if (k >= ctx->counts.n_bh || ctx->bh_all[k].core.swallow_id != og3_debug_marker64(id, 3u)) errors++;
        }
        if (errors) {
            fprintf(stderr, "Bc deep verify failed at slot %llu\n", (unsigned long long)i);
            return errors;
        }
    }
    return 0;
#endif
}
