/* =========================================================================
 * layout_C.c
 *
 * Layout C: fine common atoms + fine type-specific arrays + positional
 * mapping. This C benchmark intentionally keeps a resident type_idx[] side
 * array so same-type swaps and verification can locate the implied type slot
 * cheaply; its memory/refresh cost is reported as bookkeeping.
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

    PCore *core;
    PDyn  *dyn;
    PTime *time;
    PMeta *meta;

    GasCore  *gas_core;
    GasGrad  *gas_grad;
    GasMetal *gas_metal;
    GasSF    *gas_sf;

    StarCore *star_core;
    StarMeta *star_meta;

    BHCore  *bh_core;
    BHEnv   *bh_env;
    BHRepos *bh_repos;

    count_t *type_idx;
    count_t *sub_perm_gas;
    count_t *sub_perm_star;
    count_t *sub_perm_bh;

    void    *scratch;
    size_t   scratch_bytes;
    count_t *scratch_idx;
};

const char *layout_name(void)        { return "C"; }
const char *layout_description(void) { return "Fine common + fine type-specific + positional mapping"; }

layout_ctx_t *layout_alloc(const pcount_t * restrict c)
{
    layout_ctx_t *ctx = (layout_ctx_t *)calloc(1, sizeof(*ctx));
    if (!ctx) return NULL;
    ctx->counts = *c;
    count_t N = total_parts(c);

    ctx->core = (PCore *)og3_aligned_alloc((size_t)N * sizeof(PCore));
    ctx->dyn  = (PDyn  *)og3_aligned_alloc((size_t)N * sizeof(PDyn));
    ctx->time = (PTime *)og3_aligned_alloc((size_t)N * sizeof(PTime));
    ctx->meta = (PMeta *)og3_aligned_alloc((size_t)N * sizeof(PMeta));
    ctx->type_idx = (count_t *)og3_aligned_alloc((size_t)N * sizeof(count_t));

    if (c->n_gas) {
        ctx->gas_core  = (GasCore  *)og3_aligned_alloc((size_t)c->n_gas * sizeof(GasCore));
        ctx->gas_grad  = (GasGrad  *)og3_aligned_alloc((size_t)c->n_gas * sizeof(GasGrad));
        ctx->gas_metal = (GasMetal *)og3_aligned_alloc((size_t)c->n_gas * sizeof(GasMetal));
        ctx->gas_sf    = (GasSF    *)og3_aligned_alloc((size_t)c->n_gas * sizeof(GasSF));
        ctx->sub_perm_gas = (count_t *)og3_aligned_alloc((size_t)c->n_gas * sizeof(count_t));
    }
    if (c->n_star) {
        ctx->star_core = (StarCore *)og3_aligned_alloc((size_t)c->n_star * sizeof(StarCore));
        ctx->star_meta = (StarMeta *)og3_aligned_alloc((size_t)c->n_star * sizeof(StarMeta));
        ctx->sub_perm_star = (count_t *)og3_aligned_alloc((size_t)c->n_star * sizeof(count_t));
    }
    if (c->n_bh) {
        ctx->bh_core  = (BHCore  *)og3_aligned_alloc((size_t)c->n_bh * sizeof(BHCore));
        ctx->bh_env   = (BHEnv   *)og3_aligned_alloc((size_t)c->n_bh * sizeof(BHEnv));
        ctx->bh_repos = (BHRepos *)og3_aligned_alloc((size_t)c->n_bh * sizeof(BHRepos));
        ctx->sub_perm_bh = (count_t *)og3_aligned_alloc((size_t)c->n_bh * sizeof(count_t));
    }

    size_t max_atom = sizeof(PCore);
    if (sizeof(PDyn)     > max_atom) max_atom = sizeof(PDyn);
    if (sizeof(PTime)    > max_atom) max_atom = sizeof(PTime);
    if (sizeof(PMeta)    > max_atom) max_atom = sizeof(PMeta);
    if (sizeof(GasCore)  > max_atom) max_atom = sizeof(GasCore);
    if (sizeof(GasGrad)  > max_atom) max_atom = sizeof(GasGrad);
    if (sizeof(GasMetal) > max_atom) max_atom = sizeof(GasMetal);
    if (sizeof(GasSF)    > max_atom) max_atom = sizeof(GasSF);
    if (sizeof(StarCore) > max_atom) max_atom = sizeof(StarCore);
    if (sizeof(StarMeta) > max_atom) max_atom = sizeof(StarMeta);
    if (sizeof(BHCore)   > max_atom) max_atom = sizeof(BHCore);
    if (sizeof(BHEnv)    > max_atom) max_atom = sizeof(BHEnv);
    if (sizeof(BHRepos)  > max_atom) max_atom = sizeof(BHRepos);
    ctx->scratch_bytes = (size_t)N * max_atom;
    ctx->scratch = og3_aligned_alloc(ctx->scratch_bytes);
    ctx->scratch_idx = (count_t *)og3_aligned_alloc((size_t)N * sizeof(count_t));

    if ((N && (!ctx->core || !ctx->dyn || !ctx->time || !ctx->meta ||
               !ctx->type_idx || !ctx->scratch || !ctx->scratch_idx)) ||
        (c->n_gas && (!ctx->gas_core || !ctx->gas_grad || !ctx->gas_metal || !ctx->gas_sf || !ctx->sub_perm_gas)) ||
        (c->n_star && (!ctx->star_core || !ctx->star_meta || !ctx->sub_perm_star)) ||
        (c->n_bh && (!ctx->bh_core || !ctx->bh_env || !ctx->bh_repos || !ctx->sub_perm_bh))) {
        layout_free(ctx);
        return NULL;
    }
    return ctx;
}

void layout_free(layout_ctx_t * restrict ctx)
{
    if (!ctx) return;
    free(ctx->core); free(ctx->dyn); free(ctx->time); free(ctx->meta);
    free(ctx->gas_core); free(ctx->gas_grad); free(ctx->gas_metal); free(ctx->gas_sf);
    free(ctx->star_core); free(ctx->star_meta);
    free(ctx->bh_core); free(ctx->bh_env); free(ctx->bh_repos);
    free(ctx->type_idx); free(ctx->sub_perm_gas); free(ctx->sub_perm_star); free(ctx->sub_perm_bh);
    free(ctx->scratch); free(ctx->scratch_idx);
    free(ctx);
}

static void refresh_type_idx(layout_ctx_t * restrict ctx)
{
    count_t N = total_parts(&ctx->counts);
    count_t counters[NTYPES] = {0};
    for (count_t p = 0; p < N; ++p) {
        uint8_t t = ctx->core[p].type;
        ctx->type_idx[p] = counters[t]++;
    }
    assert(counters[PT_GAS] == ctx->counts.n_gas);
    assert(counters[PT_STAR] == ctx->counts.n_star);
    assert(counters[PT_BH] == ctx->counts.n_bh);
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

        switch (t) {
        case PT_GAS:
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
            memset(&ctx->star_core[is], 0, sizeof(StarCore));
            memset(&ctx->star_meta[is], 0, sizeof(StarMeta));
            ctx->star_core[is].stellar_age = (real_t)(p * 1e-5);
#ifdef DEBUG
            ctx->star_core[is].pid = og3_debug_marker32(ids[p], 2u);
#endif
            ++is;
            break;
        case PT_BH:
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
            break;
        }
    }
    assert(ig == ctx->counts.n_gas);
    assert(is == ctx->counts.n_star);
    assert(ib == ctx->counts.n_bh);
    refresh_type_idx(ctx);
}

void layout_reshuffle_full(layout_ctx_t * restrict ctx, const count_t * restrict perm)
{
    count_t N = total_parts(&ctx->counts);

    count_t sub_counters[NTYPES] = {0};
    for (count_t i = 0; i < N; ++i) {
        count_t src = perm[i];
        uint8_t t = ctx->core[src].type;
        count_t k = ctx->type_idx[src];
        switch (t) {
        case PT_GAS:  ctx->sub_perm_gas [sub_counters[t]++] = k; break;
        case PT_STAR: ctx->sub_perm_star[sub_counters[t]++] = k; break;
        case PT_BH:   ctx->sub_perm_bh  [sub_counters[t]++] = k; break;
        default: sub_counters[t]++; break;
        }
    }
    assert(sub_counters[PT_GAS] == ctx->counts.n_gas);
    assert(sub_counters[PT_STAR] == ctx->counts.n_star);
    assert(sub_counters[PT_BH] == ctx->counts.n_bh);

    og3_apply_perm_bytes(ctx->core, perm, N, sizeof(PCore), ctx->scratch, ctx->scratch_idx);
    og3_apply_perm_bytes(ctx->dyn,  perm, N, sizeof(PDyn),  ctx->scratch, ctx->scratch_idx);
    og3_apply_perm_bytes(ctx->time, perm, N, sizeof(PTime), ctx->scratch, ctx->scratch_idx);
    og3_apply_perm_bytes(ctx->meta, perm, N, sizeof(PMeta), ctx->scratch, ctx->scratch_idx);

    if (ctx->counts.n_gas) {
        count_t Ng = ctx->counts.n_gas;
        og3_apply_perm_bytes(ctx->gas_core,  ctx->sub_perm_gas, Ng, sizeof(GasCore),  ctx->scratch, ctx->scratch_idx);
        og3_apply_perm_bytes(ctx->gas_grad,  ctx->sub_perm_gas, Ng, sizeof(GasGrad),  ctx->scratch, ctx->scratch_idx);
        og3_apply_perm_bytes(ctx->gas_metal, ctx->sub_perm_gas, Ng, sizeof(GasMetal), ctx->scratch, ctx->scratch_idx);
        og3_apply_perm_bytes(ctx->gas_sf,    ctx->sub_perm_gas, Ng, sizeof(GasSF),    ctx->scratch, ctx->scratch_idx);
    }
    if (ctx->counts.n_star) {
        count_t Ns = ctx->counts.n_star;
        og3_apply_perm_bytes(ctx->star_core, ctx->sub_perm_star, Ns, sizeof(StarCore), ctx->scratch, ctx->scratch_idx);
        og3_apply_perm_bytes(ctx->star_meta, ctx->sub_perm_star, Ns, sizeof(StarMeta), ctx->scratch, ctx->scratch_idx);
    }
    if (ctx->counts.n_bh) {
        count_t Nb = ctx->counts.n_bh;
        og3_apply_perm_bytes(ctx->bh_core,  ctx->sub_perm_bh, Nb, sizeof(BHCore),  ctx->scratch, ctx->scratch_idx);
        og3_apply_perm_bytes(ctx->bh_env,   ctx->sub_perm_bh, Nb, sizeof(BHEnv),   ctx->scratch, ctx->scratch_idx);
        og3_apply_perm_bytes(ctx->bh_repos, ctx->sub_perm_bh, Nb, sizeof(BHRepos), ctx->scratch, ctx->scratch_idx);
    }

    refresh_type_idx(ctx);
}

#define SWAP_STRUCT(a, b, S) do {    \
    S tmp;                           \
    memcpy(&tmp, &(a), sizeof(S));   \
    memcpy(&(a), &(b), sizeof(S));   \
    memcpy(&(b), &tmp, sizeof(S));   \
} while (0)

void layout_swap_same_type(layout_ctx_t * restrict ctx, count_t i, count_t j)
{
    uint8_t t = ctx->core[i].type;
    count_t ti = ctx->type_idx[i];
    count_t tj = ctx->type_idx[j];

    SWAP_STRUCT(ctx->core[i], ctx->core[j], PCore);
    SWAP_STRUCT(ctx->dyn[i],  ctx->dyn[j],  PDyn);
    SWAP_STRUCT(ctx->time[i], ctx->time[j], PTime);
    SWAP_STRUCT(ctx->meta[i], ctx->meta[j], PMeta);

    switch (t) {
    case PT_GAS:
        SWAP_STRUCT(ctx->gas_core[ti],  ctx->gas_core[tj],  GasCore);
        SWAP_STRUCT(ctx->gas_grad[ti],  ctx->gas_grad[tj],  GasGrad);
        SWAP_STRUCT(ctx->gas_metal[ti], ctx->gas_metal[tj], GasMetal);
        SWAP_STRUCT(ctx->gas_sf[ti],    ctx->gas_sf[tj],    GasSF);
        break;
    case PT_STAR:
        SWAP_STRUCT(ctx->star_core[ti], ctx->star_core[tj], StarCore);
        SWAP_STRUCT(ctx->star_meta[ti], ctx->star_meta[tj], StarMeta);
        break;
    case PT_BH:
        SWAP_STRUCT(ctx->bh_core[ti],  ctx->bh_core[tj],  BHCore);
        SWAP_STRUCT(ctx->bh_env[ti],   ctx->bh_env[tj],   BHEnv);
        SWAP_STRUCT(ctx->bh_repos[ti], ctx->bh_repos[tj], BHRepos);
        break;
    default:
        break;
    }
}

size_t layout_reshuffle_bytes(const layout_ctx_t * restrict ctx)
{
    count_t N = total_parts(&ctx->counts);
    size_t common = (size_t)N * (sizeof(PCore) + sizeof(PDyn) + sizeof(PTime) + sizeof(PMeta));
    size_t gas = (size_t)ctx->counts.n_gas * (sizeof(GasCore) + sizeof(GasGrad) + sizeof(GasMetal) + sizeof(GasSF));
    size_t star = (size_t)ctx->counts.n_star * (sizeof(StarCore) + sizeof(StarMeta));
    size_t bh = (size_t)ctx->counts.n_bh * (sizeof(BHCore) + sizeof(BHEnv) + sizeof(BHRepos));
    size_t type_idx_refresh = (size_t)N * sizeof(count_t);
    return og3_perm_dense_copy_bytes(1, common + gas + star + bh) + type_idx_refresh;
}

int layout_reshuffle_streams(const layout_ctx_t * restrict ctx)
{
    int n = 4;
    if (ctx->counts.n_gas) n += 4;
    if (ctx->counts.n_star) n += 2;
    if (ctx->counts.n_bh) n += 3;
    return n;
}

size_t layout_swap_bytes(const layout_ctx_t * restrict ctx, uint8_t t)
{
    (void)ctx;
    size_t common = sizeof(PCore) + sizeof(PDyn) + sizeof(PTime) + sizeof(PMeta);
    size_t type = 0;
    switch (t) {
    case PT_GAS:  type = sizeof(GasCore) + sizeof(GasGrad) + sizeof(GasMetal) + sizeof(GasSF); break;
    case PT_STAR: type = sizeof(StarCore) + sizeof(StarMeta); break;
    case PT_BH:   type = sizeof(BHCore) + sizeof(BHEnv) + sizeof(BHRepos); break;
    default: break;
    }
    return 6ull * (common + type);
}

int layout_swap_streams(const layout_ctx_t * restrict ctx, uint8_t t)
{
    (void)ctx;
    int n = 4;
    switch (t) {
    case PT_GAS: n += 4; break;
    case PT_STAR: n += 2; break;
    case PT_BH: n += 3; break;
    default: break;
    }
    return n;
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
        count_t k = ctx->type_idx[i];
        if (t == PT_GAS) {
            if (k >= ctx->counts.n_gas || ctx->gas_sf[k].swallow_id != og3_debug_marker64(id, 1u)) errors++;
        } else if (t == PT_STAR) {
            if (k >= ctx->counts.n_star || ctx->star_core[k].pid != og3_debug_marker32(id, 2u)) errors++;
        } else if (t == PT_BH) {
            if (k >= ctx->counts.n_bh || ctx->bh_core[k].swallow_id != og3_debug_marker64(id, 3u)) errors++;
        }
        if (errors) {
            fprintf(stderr, "C deep verify failed at slot %llu\n", (unsigned long long)i);
            return errors;
        }
    }
    return 0;
#endif
}
