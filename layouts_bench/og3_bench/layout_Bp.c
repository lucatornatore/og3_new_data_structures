/* =========================================================================
 * layout_Bp.c
 *
 * Layout B': fine common atoms + fine type-specific arrays + explicit
 * cross-index. Global PH reshuffle permutes only the common arrays; type
 * arrays stay put and are reached through PLinkage.type_idx.
 * ========================================================================= */

#include "layout_api.h"
#include "fields.h"
#include "bench_common.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct layout_ctx {
    pcount_t counts;

    PCore    *core;
    PDyn     *dyn;
    PTime    *time;
    PMeta    *meta;
    PLinkage *linkage;

    GasCore  *gas_core;
    GasGrad  *gas_grad;
    GasMetal *gas_metal;
    GasSF    *gas_sf;

    StarCore *star_core;
    StarMeta *star_meta;

    BHCore   *bh_core;
    BHEnv    *bh_env;
    BHRepos  *bh_repos;

    void     *scratch;
    size_t    scratch_bytes;
    count_t  *scratch_idx;
};

const char *layout_name(void)        { return "B'"; }
const char *layout_description(void) { return "Fine common + fine type-specific + cross-index"; }

layout_ctx_t *layout_alloc(const pcount_t * restrict c)
{
    layout_ctx_t *ctx = (layout_ctx_t *)calloc(1, sizeof(*ctx));
    if (!ctx) return NULL;
    ctx->counts = *c;
    count_t N = total_parts(c);

    ctx->core    = (PCore    *)og3_aligned_alloc((size_t)N * sizeof(PCore));
    ctx->dyn     = (PDyn     *)og3_aligned_alloc((size_t)N * sizeof(PDyn));
    ctx->time    = (PTime    *)og3_aligned_alloc((size_t)N * sizeof(PTime));
    ctx->meta    = (PMeta    *)og3_aligned_alloc((size_t)N * sizeof(PMeta));
    ctx->linkage = (PLinkage *)og3_aligned_alloc((size_t)N * sizeof(PLinkage));

    if (c->n_gas) {
        ctx->gas_core  = (GasCore  *)og3_aligned_alloc((size_t)c->n_gas * sizeof(GasCore));
        ctx->gas_grad  = (GasGrad  *)og3_aligned_alloc((size_t)c->n_gas * sizeof(GasGrad));
        ctx->gas_metal = (GasMetal *)og3_aligned_alloc((size_t)c->n_gas * sizeof(GasMetal));
        ctx->gas_sf    = (GasSF    *)og3_aligned_alloc((size_t)c->n_gas * sizeof(GasSF));
    }
    if (c->n_star) {
        ctx->star_core = (StarCore *)og3_aligned_alloc((size_t)c->n_star * sizeof(StarCore));
        ctx->star_meta = (StarMeta *)og3_aligned_alloc((size_t)c->n_star * sizeof(StarMeta));
    }
    if (c->n_bh) {
        ctx->bh_core  = (BHCore  *)og3_aligned_alloc((size_t)c->n_bh * sizeof(BHCore));
        ctx->bh_env   = (BHEnv   *)og3_aligned_alloc((size_t)c->n_bh * sizeof(BHEnv));
        ctx->bh_repos = (BHRepos *)og3_aligned_alloc((size_t)c->n_bh * sizeof(BHRepos));
    }

    size_t max_atom = sizeof(PCore);
    if (sizeof(PDyn)     > max_atom) max_atom = sizeof(PDyn);
    if (sizeof(PTime)    > max_atom) max_atom = sizeof(PTime);
    if (sizeof(PMeta)    > max_atom) max_atom = sizeof(PMeta);
    if (sizeof(PLinkage) > max_atom) max_atom = sizeof(PLinkage);
    ctx->scratch_bytes = (size_t)N * max_atom;
    ctx->scratch = og3_aligned_alloc(ctx->scratch_bytes);
    ctx->scratch_idx = (count_t *)og3_aligned_alloc((size_t)N * sizeof(count_t));

    if ((N && (!ctx->core || !ctx->dyn || !ctx->time || !ctx->meta || !ctx->linkage ||
               !ctx->scratch || !ctx->scratch_idx)) ||
        (c->n_gas && (!ctx->gas_core || !ctx->gas_grad || !ctx->gas_metal || !ctx->gas_sf)) ||
        (c->n_star && (!ctx->star_core || !ctx->star_meta)) ||
        (c->n_bh && (!ctx->bh_core || !ctx->bh_env || !ctx->bh_repos))) {
        layout_free(ctx);
        return NULL;
    }
    return ctx;
}

void layout_free(layout_ctx_t * restrict ctx)
{
    if (!ctx) return;
    free(ctx->core); free(ctx->dyn); free(ctx->time); free(ctx->meta); free(ctx->linkage);
    free(ctx->gas_core); free(ctx->gas_grad); free(ctx->gas_metal); free(ctx->gas_sf);
    free(ctx->star_core); free(ctx->star_meta);
    free(ctx->bh_core); free(ctx->bh_env); free(ctx->bh_repos);
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

        memset(&ctx->core[p], 0, sizeof(PCore));
        ctx->core[p].key = keys[p];
        ctx->core[p].type = t;
        ctx->core[p].pos[0] = (pos_t)(p * 1e-6);

        memset(&ctx->dyn[p], 0, sizeof(PDyn));
        ctx->dyn[p].vel[0] = (real_t)(p * 1e-3);

        memset(&ctx->time[p], 0, sizeof(PTime));
        ctx->time[p].ti_current = (time_int_t)p;

        memset(&ctx->meta[p], 0, sizeof(PMeta));
        ctx->meta[p].id = ids[p];

        memset(&ctx->linkage[p], 0, sizeof(PLinkage));
        switch (t) {
        case PT_GAS:
            ctx->linkage[p].type_idx = ig;
            memset(&ctx->gas_core[ig], 0, sizeof(GasCore));
            memset(&ctx->gas_grad[ig], 0, sizeof(GasGrad));
            memset(&ctx->gas_metal[ig], 0, sizeof(GasMetal));
            memset(&ctx->gas_sf[ig], 0, sizeof(GasSF));
            ctx->gas_core[ig].density = (real_t)(p * 1e-4);
            ctx->gas_metal[ig].metals[0] = (real_t)(p % 1000);
#ifdef DEBUG
            ctx->gas_sf[ig].swallow_id = og3_debug_marker64(ids[p], 1u);
#endif
            ++ig;
            break;
        case PT_STAR:
            ctx->linkage[p].type_idx = is;
            memset(&ctx->star_core[is], 0, sizeof(StarCore));
            memset(&ctx->star_meta[is], 0, sizeof(StarMeta));
            ctx->star_core[is].stellar_age = (real_t)(p * 1e-5);
#ifdef DEBUG
            ctx->star_core[is].pid = og3_debug_marker32(ids[p], 2u);
#endif
            ++is;
            break;
        case PT_BH:
            ctx->linkage[p].type_idx = ib;
            memset(&ctx->bh_core[ib], 0, sizeof(BHCore));
            memset(&ctx->bh_env[ib], 0, sizeof(BHEnv));
            memset(&ctx->bh_repos[ib], 0, sizeof(BHRepos));
            ctx->bh_core[ib].bh_mass = (double)(p * 1e-6);
#ifdef DEBUG
            ctx->bh_core[ib].swallow_id = og3_debug_marker64(ids[p], 3u);
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
    og3_apply_perm_bytes(ctx->core,    perm, N, sizeof(PCore),    ctx->scratch, ctx->scratch_idx);
    og3_apply_perm_bytes(ctx->dyn,     perm, N, sizeof(PDyn),     ctx->scratch, ctx->scratch_idx);
    og3_apply_perm_bytes(ctx->time,    perm, N, sizeof(PTime),    ctx->scratch, ctx->scratch_idx);
    og3_apply_perm_bytes(ctx->meta,    perm, N, sizeof(PMeta),    ctx->scratch, ctx->scratch_idx);
    og3_apply_perm_bytes(ctx->linkage, perm, N, sizeof(PLinkage), ctx->scratch, ctx->scratch_idx);
}

#define SWAP_STRUCT(a, b, S) do { \
    S tmp;                       \
    memcpy(&tmp, &(a), sizeof(S)); \
    memcpy(&(a), &(b), sizeof(S)); \
    memcpy(&(b), &tmp, sizeof(S)); \
} while (0)

void layout_swap_same_type(layout_ctx_t * restrict ctx, count_t i, count_t j)
{
    SWAP_STRUCT(ctx->core[i],    ctx->core[j],    PCore);
    SWAP_STRUCT(ctx->dyn[i],     ctx->dyn[j],     PDyn);
    SWAP_STRUCT(ctx->time[i],    ctx->time[j],    PTime);
    SWAP_STRUCT(ctx->meta[i],    ctx->meta[j],    PMeta);
    SWAP_STRUCT(ctx->linkage[i], ctx->linkage[j], PLinkage);
}

size_t layout_reshuffle_bytes(const layout_ctx_t * restrict ctx)
{
    count_t N = total_parts(&ctx->counts);
    size_t per = sizeof(PCore) + sizeof(PDyn) + sizeof(PTime) + sizeof(PMeta) + sizeof(PLinkage);
    return og3_perm_dense_copy_bytes(N, per);
}

int layout_reshuffle_streams(const layout_ctx_t * restrict ctx)
{
    (void)ctx;
    return 5;
}

size_t layout_swap_bytes(const layout_ctx_t * restrict ctx, uint8_t t)
{
    (void)ctx; (void)t;
    return 6ull * (sizeof(PCore) + sizeof(PDyn) + sizeof(PTime) + sizeof(PMeta) + sizeof(PLinkage));
}

int layout_swap_streams(const layout_ctx_t * restrict ctx, uint8_t t)
{
    (void)ctx; (void)t;
    return 5;
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
ppid_t layout_get_id(const layout_ctx_t * restrict ctx, count_t i) { return ctx->meta[i].id; }

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
        ppid_t id = ctx->meta[i].id;
        count_t k = ctx->linkage[i].type_idx;
        if (t == PT_GAS) {
            if (k >= ctx->counts.n_gas || ctx->gas_sf[k].swallow_id != og3_debug_marker64(id, 1u)) errors++;
        } else if (t == PT_STAR) {
            if (k >= ctx->counts.n_star || ctx->star_core[k].pid != og3_debug_marker32(id, 2u)) errors++;
        } else if (t == PT_BH) {
            if (k >= ctx->counts.n_bh || ctx->bh_core[k].swallow_id != og3_debug_marker64(id, 3u)) errors++;
        }
        if (errors) {
            fprintf(stderr, "B' deep verify failed at slot %llu\n", (unsigned long long)i);
            return errors;
        }
    }
    return 0;
#endif
}
