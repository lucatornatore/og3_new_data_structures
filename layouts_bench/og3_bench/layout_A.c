/* =========================================================================
 * layout_A.c
 *
 * LAYOUT A — Monolithic per-type AoS.
 *
 * One fat struct per particle type. No cross-indexing: common fields
 * (pos, vel, key, ...) are DUPLICATED in each per-type struct, because
 * types never change under a pure reshuffle (migration and mutation are
 * the only operations that change them, and neither is in scope here).
 *
 * Arrays:
 *   gas_[]   : N_gas  * GasFullA    (pos..sound_speed..metals..sfr)
 *   dm_[]    : N_dm   * DMFullA     (pos..vel..grav_acc..id — common only)
 *   star_[]  : N_star * StarFullA   (common + star fields)
 *   bh_[]    : N_bh   * BHFullA     (common + bh fields)
 *
 * "Common" fields are those that apply to every particle type:
 *   pos, mass, key, type, flags, time_bin, leaf_idx,
 *   vel, grav_acc, old_acc, num_ngb, true_ngb,
 *   ti_begin, ti_current, dt_step, grav_cost,
 *   id
 * Sum: ~128 B before padding.
 *
 * Benchmarks exercised on Layout A have a unique character:
 *   - The benchmark builds ONE global logical PH-key permutation. A cannot
 *     physically apply that permutation to one common array because it owns
 *     four per-type fat arrays. Therefore layout_reshuffle_full() first
 *     permutes the logical SlotA[] directory and derives four per-type
 *     sub-permutations from the global order.
 *   - Swap is a single fat-struct swap inside one per-type array, after
 *     translating global SlotA indices to type-local indices.
 *   - Bytes moved per reshuffle is large (fat structs) but streams are
 *     few (four independent arrays plus the SlotA[] directory copy).
 * ========================================================================= */

#include "layout_api.h"
#include "bench_common.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>

/* ------------------------------------------------------------------------ */
/* Fat per-type struct layouts                                              */
/*                                                                          */
/* Each per-type struct inlines the same "common" prefix (same field layout */
/* and offsets) so that fill_common_prefix can memcpy into any of them.     */
/* This prefix matches the atom decomposition used by Bp/Bc/C:              */
/*   PCore + PDyn + PTime + PMeta  =  atoms of B'/Bc/C                      */
/* Then the type-specific suffix inlines the type-specific atoms.           */
/* ------------------------------------------------------------------------ */

/* Shared common prefix — identical leading layout in all four fat structs.
 * Natural alignment walk (verified below):
 *   pos[3]     0..24
 *   mass       24..28
 *   key        32..40    (pad 4 before, 8-aligned)
 *   type,flags,time_bin,leaf_idx  40..48
 *   vel[3]     48..60
 *   dp[3]      60..72
 *   grav_acc[3] 72..96   (8-aligned OK)
 *   grav_pm[3] 96..120
 *   old_acc    120..124
 *   num_ngb    124..128
 *   true_ngb   128..132
 *   ti_begin   136..144  (pad 4 before, 8-aligned)
 *   ti_current 144..152
 *   dt_step    152..156
 *   grav_cost[6] 156..180
 *   id         184..192  (pad 4 before, 8-aligned)
 *   hsml_index 192..196
 *   flags_ext  196..200
 * common raw = 200 B.
 */
#define A_COMMON_FIELDS                                                     \
    pos_t      pos[3];                                                      \
    real_t     mass;                                                        \
    pkey_t     key;                                                         \
    uint8_t    type;                                                        \
    uint8_t    flags;                                                       \
    int16_t    time_bin;                                                    \
    int32_t    leaf_idx;                                                    \
    real_t     vel[3];                                                      \
    real_t     dp[3];                                                       \
    acc_t      grav_acc[3];                                                 \
    acc_t      grav_pm[3];                                                  \
    real_t     old_acc;                                                     \
    real_t     num_ngb;                                                     \
    int32_t    true_ngb;                                                    \
    time_int_t ti_begin;                                                    \
    time_int_t ti_current;                                                  \
    int32_t    dt_step;                                                     \
    real_t     grav_cost[6];                                                \
    ppid_t     id;                                                          \
    uint32_t   hsml_index;                                                  \
    uint32_t   flags_ext;

/* Gas: adds GasCore+GasGrad+GasMetal+GasSF fields inline. */
typedef struct OG3_ATOM {
    A_COMMON_FIELDS
    /* GasCore (60 B raw) */
    real_t   density;
    real_t   density_old;
    real_t   pressure;
    real_t   entropy;
    real_t   entropy_pred;
    real_t   dt_entropy;
    acc_t    dhsml_factor;
    real_t   alpha;
    real_t   f_balsara;
    real_t   temperature;
    real_t   x_cold_cloud;
    real_t   elec;
    real_t   injected_bh_energy;
    real_t   max_signal_vel;
    /* GasGrad (64 B raw) */
    real_t   vel_pred[3];
    acc_t    hydro_acc[3];
    real_t   div_v;
    real_t   curl_v;
    acc_t    egy_step;
    double   mstar;
    /* GasMetal (64 B raw) */
    real_t   metals[16];
    /* GasSF (24 B raw) */
    ppid_t   swallow_id;
    real_t   sfr;
    real_t   delay_time;
    real_t   mass_res;
    real_t   egy_res;
    /* raw ~408; aligned(64) => 448 */
} GasFullA;

/* DM: no type-specific fields. */
typedef struct OG3_ATOM {
    A_COMMON_FIELDS
    /* raw 200; aligned(64) => 256 */
} DMFullA;

/* Star: adds StarCore+StarMeta fields inline. */
typedef struct OG3_ATOM {
    A_COMMON_FIELDS
    /* StarCore (28 B raw) — weight is double, must be 8-aligned; offset 200 OK */
    double   star_weight;
    real_t   last_chem_time;
    real_t   i_mass;
    uint32_t pid_star;
    real_t   stellar_age;
    int32_t  chem_time_bin;
    /* StarMeta (64 B raw) */
    real_t   star_metals[16];
    /* raw ~292; aligned(64) => 320 */
} StarFullA;

/* BH: adds BHCore+BHEnv+BHRepos fields inline. */
typedef struct OG3_ATOM {
    A_COMMON_FIELDS
    /* BHCore (48 B raw) */
    uint32_t pid_bh;
    real_t   stellar_age_bh;
    ppid_t   swallow_id_bh;
    real_t   swallow_pot;
    int32_t  count_progs;
    double   bh_mass;
    real_t   bh_mdot;
    int32_t  time_bin_gas_ngb;
    real_t   bh_density;
    real_t   bh_entropy;
    /* BHEnv (52 B raw) */
    real_t   surr_gas_vel[3];
    real_t   cold_density;
    real_t   cold_entropy;
    real_t   surr_cold_vel[3];
    real_t   hot_density;
    real_t   hot_entropy;
    real_t   surr_hot_vel[3];
    /* BHRepos (116 B raw) */
    real_t   accreted_mass;
    real_t   accreted_bh_mass;
    real_t   accreted_mom[3];
    int32_t  flag_merged;
    double   final_spin[3];
    real_t   surr_vel[3];
    real_t   surr_density;
    real_t   sigma;
    real_t   bmax;
    real_t   tot_fb_eff;
    real_t   mean_hsml;
    real_t   mean_rho;
    double   min_pot_pos[3];
    real_t   min_pot;
    /* raw ~420; aligned(64) => 448 */
} BHFullA;

typedef struct {
    uint8_t type;
    count_t type_idx;
} SlotA;

/* ------------------------------------------------------------------------ */
/* Layout A context                                                         */
/* ------------------------------------------------------------------------ */

struct layout_ctx {
    pcount_t    counts;
    GasFullA   *gas;
    DMFullA    *dm;
    StarFullA  *star;
    BHFullA    *bh;

    /* Logical global PH order. Slot i tells which per-type fat array owns
     * the particle that occupies global slot i. */
    SlotA      *slots;
    SlotA      *slot_scratch;

    count_t    *scratch;      /* sub-permutation workspace */
    count_t    *scratch_idx;  /* changed-index workspace for sparse apply */
    void       *scratch_data; /* scratch for struct-sized reshuffle */
    size_t      scratch_data_bytes;
};

/* ------------------------------------------------------------------------ */
/* API                                                                      */
/* ------------------------------------------------------------------------ */

const char *layout_name(void)        { return "A"; }
const char *layout_description(void) { return "Monolithic per-type AoS (fat struct, no cross-index)"; }

layout_ctx_t *layout_alloc(const pcount_t * restrict c)
{
    layout_ctx_t *ctx = (layout_ctx_t *)calloc(1, sizeof(layout_ctx_t));
    if (!ctx) return NULL;

    ctx->counts = *c;

    if (c->n_gas)  ctx->gas  = (GasFullA  *)og3_aligned_alloc(c->n_gas  * sizeof(GasFullA));
    if (c->n_dm)   ctx->dm   = (DMFullA   *)og3_aligned_alloc(c->n_dm   * sizeof(DMFullA));
    if (c->n_star) ctx->star = (StarFullA *)og3_aligned_alloc(c->n_star * sizeof(StarFullA));
    if (c->n_bh)   ctx->bh   = (BHFullA   *)og3_aligned_alloc(c->n_bh   * sizeof(BHFullA));

    /* Scratch for reshuffle: max struct size * max per-type count. */
    count_t max_n   = c->n_gas;
    size_t  max_sz  = sizeof(GasFullA);
    if (c->n_dm   > max_n) max_n = c->n_dm;
    if (c->n_star > max_n) max_n = c->n_star;
    if (c->n_bh   > max_n) max_n = c->n_bh;
    if (sizeof(DMFullA)   > max_sz) max_sz = sizeof(DMFullA);
    if (sizeof(StarFullA) > max_sz) max_sz = sizeof(StarFullA);
    if (sizeof(BHFullA)   > max_sz) max_sz = sizeof(BHFullA);

    count_t N = total_parts(c);
    ctx->scratch_data       = og3_aligned_alloc((size_t)max_n * max_sz);
    ctx->scratch_data_bytes = (size_t)max_n * max_sz;
    ctx->scratch            = (count_t *)og3_aligned_alloc((size_t)N * sizeof(count_t));
    ctx->scratch_idx        = (count_t *)og3_aligned_alloc((size_t)N * sizeof(count_t));
    ctx->slots              = (SlotA *)og3_aligned_alloc((size_t)N * sizeof(SlotA));
    ctx->slot_scratch       = (SlotA *)og3_aligned_alloc((size_t)N * sizeof(SlotA));

    if ((!ctx->gas && c->n_gas) || (!ctx->dm && c->n_dm) ||
        (!ctx->star && c->n_star) || (!ctx->bh && c->n_bh) ||
        (N && (!ctx->scratch_data || !ctx->scratch || !ctx->scratch_idx || !ctx->slots || !ctx->slot_scratch))) {
        layout_free(ctx);
        return NULL;
    }
    return ctx;
}

void layout_free(layout_ctx_t * restrict ctx)
{
    if (!ctx) return;
    free(ctx->gas);
    free(ctx->dm);
    free(ctx->star);
    free(ctx->bh);
    free(ctx->scratch);
    free(ctx->scratch_idx);
    free(ctx->scratch_data);
    free(ctx->slots);
    free(ctx->slot_scratch);
    free(ctx);
}

/* ------------------------------------------------------------------------
 * fill: walks the caller-supplied arrays and distributes particles into
 * the four per-type arrays in the order they are encountered. The
 * per-type index assigned is simply the running count.
 * ------------------------------------------------------------------------ */

/* helper: write common fields to any of the four structs. Because the
 * four structs have IDENTICAL leading layout for the common region (this
 * is enforced by Static-assertion of offset equality below), we can
 * memcpy the prefix. */

_Static_assert(offsetof(DMFullA,   id) == offsetof(GasFullA, id),
               "common prefix must be positionally identical across A structs");
_Static_assert(offsetof(StarFullA, id) == offsetof(GasFullA, id),
               "common prefix must be positionally identical across A structs");
_Static_assert(offsetof(BHFullA,   id) == offsetof(GasFullA, id),
               "common prefix must be positionally identical across A structs");

static void fill_common_prefix(void *dst, pkey_t key, uint8_t type, ppid_t id,
                               count_t p)
{
    /* We only set the fields a benchmark actually exercises: key, type,
     * id. Everything else is zero-initialised. */
    GasFullA tmp = {0};
    tmp.key       = key;
    tmp.type      = type;
    tmp.id        = id;
    tmp.pos[0]    = (pos_t)(p * 1e-6);
    tmp.vel[0]    = (real_t)(p * 1e-3);
    /* Common region = everything before the first gas-specific field.
     * offsetof(GasFullA, density) marks that boundary. */
    memcpy(dst, &tmp, offsetof(GasFullA, density));
}

void layout_fill(layout_ctx_t *restrict ctx,
                 const pkey_t *restrict keys, const uint8_t *restrict types, const ppid_t *restrict ids)
{
    count_t N  = total_parts(&ctx->counts);
    count_t ig = 0, idm = 0, is = 0, ib = 0;

    for (count_t p = 0; p < N; ++p) {
        uint8_t t = types[p];
        pkey_t  k = keys[p];
        ppid_t  i = ids[p];

        switch (t) {
            case PT_GAS:
                ctx->slots[p] = (SlotA){ .type = t, .type_idx = ig };
                fill_common_prefix(&ctx->gas[ig], k, t, i, p);
                ctx->gas[ig].density   = (real_t)(p * 1e-4);
                ctx->gas[ig].metals[0] = (real_t)(p % 1000);
#ifdef DEBUG
                ctx->gas[ig].swallow_id = og3_debug_marker64(i, 1u);
#endif
                ++ig;
                break;
            case PT_DM1: case PT_DM2: case PT_DM3:
                ctx->slots[p] = (SlotA){ .type = t, .type_idx = idm };
                fill_common_prefix(&ctx->dm[idm], k, t, i, p);
                ++idm;
                break;
            case PT_STAR:
                ctx->slots[p] = (SlotA){ .type = t, .type_idx = is };
                fill_common_prefix(&ctx->star[is], k, t, i, p);
                ctx->star[is].stellar_age = (real_t)(p * 1e-5);
#ifdef DEBUG
                ctx->star[is].pid_star = og3_debug_marker32(i, 2u);
#endif
                ++is;
                break;
            case PT_BH:
                ctx->slots[p] = (SlotA){ .type = t, .type_idx = ib };
                fill_common_prefix(&ctx->bh[ib], k, t, i, p);
                ctx->bh[ib].bh_mass = (double)(p * 1e-6);
#ifdef DEBUG
                ctx->bh[ib].swallow_id_bh = og3_debug_marker64(i, 3u);
#endif
                ++ib;
                break;
            default: break;
        }
    }
    assert(ig  == ctx->counts.n_gas);
    assert(idm == ctx->counts.n_dm);
    assert(is  == ctx->counts.n_star);
    assert(ib  == ctx->counts.n_bh);
}

/* ------------------------------------------------------------------------
 * reshuffle_full under A: apply one global logical permutation by
 * translating it into four per-type sub-permutations and applying each
 * sub-permutation to the corresponding fat array.
 *
 * The global perm has length N; perm[i] = j means "new slot i gets
 * what was at old slot j". For A there is no global slot ordering,
 * so we interpret perm as a specification over a hypothetical common
 * array of length N and derive per-type sub-perms by walking it in
 * order: the k-th time we see a gas particle in perm gives the new
 * position of that gas slot.
 * ------------------------------------------------------------------------ */

/* The byte-level permutation helper lives in bench_common.c. */

/* Derive per-type sub-permutations by walking perm and, for each global
 * slot, looking up what type the particle at the (OLD) slot perm[i] has.
 *
 * The problem: A does not store particle payloads in global-slot order.
 * It stores payloads in per-type arrays and stores the logical global
 * order in SlotA[]. Therefore we inspect ctx->slots[perm[i]] to ask
 * which old per-type payload should become the i-th global particle,
 * then append that payload's old type_idx to the relevant sub-perm.
 * At the same time, slot_scratch[i] is rewritten to point at the new
 * type-local position that particle will have after the type array is
 * reshuffled.
 * ------------------------------------------------------------------------ */
void layout_reshuffle_full(layout_ctx_t * restrict ctx, const count_t * restrict perm)
{
    const pcount_t *c = &ctx->counts;
    count_t N = total_parts(c);

    count_t *sub_gas  = ctx->scratch;
    count_t *sub_dm   = sub_gas  + c->n_gas;
    count_t *sub_star = sub_dm   + c->n_dm;
    count_t *sub_bh   = sub_star + c->n_star;

    count_t kg = 0, kdm = 0, ks = 0, kb = 0;
    for (count_t i = 0; i < N; ++i) {
        SlotA old = ctx->slots[perm[i]];
        switch (old.type) {
        case PT_GAS:
            sub_gas[kg] = old.type_idx;
            ctx->slot_scratch[i] = (SlotA){ .type = old.type, .type_idx = kg++ };
            break;
        case PT_DM1: case PT_DM2: case PT_DM3:
            sub_dm[kdm] = old.type_idx;
            ctx->slot_scratch[i] = (SlotA){ .type = old.type, .type_idx = kdm++ };
            break;
        case PT_STAR:
            sub_star[ks] = old.type_idx;
            ctx->slot_scratch[i] = (SlotA){ .type = old.type, .type_idx = ks++ };
            break;
        case PT_BH:
            sub_bh[kb] = old.type_idx;
            ctx->slot_scratch[i] = (SlotA){ .type = old.type, .type_idx = kb++ };
            break;
        default:
            break;
        }
    }
    assert(kg == c->n_gas);
    assert(kdm == c->n_dm);
    assert(ks == c->n_star);
    assert(kb == c->n_bh);

    if (c->n_gas)
        og3_apply_perm_bytes(ctx->gas, sub_gas, c->n_gas, sizeof(GasFullA), ctx->scratch_data, ctx->scratch_idx);
    if (c->n_dm)
        og3_apply_perm_bytes(ctx->dm, sub_dm, c->n_dm, sizeof(DMFullA), ctx->scratch_data, ctx->scratch_idx);
    if (c->n_star)
        og3_apply_perm_bytes(ctx->star, sub_star, c->n_star, sizeof(StarFullA), ctx->scratch_data, ctx->scratch_idx);
    if (c->n_bh)
        og3_apply_perm_bytes(ctx->bh, sub_bh, c->n_bh, sizeof(BHFullA), ctx->scratch_data, ctx->scratch_idx);

    memcpy(ctx->slots, ctx->slot_scratch, (size_t)N * sizeof(SlotA));
}

/* ------------------------------------------------------------------------
 * swap_same_type
 *
 * For A, same-type swap is a single-array swap of the fat struct. The
 * caller's (i, j) are global slot indices (as in the synthesised global
 * ordering used by reshuffle); we translate to per-type sub-indices.
 * ------------------------------------------------------------------------ */

#define SWAP_BYTES(p1, p2, n, tmpbuf) do {      \
    memcpy((tmpbuf), (p1),    (n));             \
    memcpy((p1),     (p2),    (n));             \
    memcpy((p2),     (tmpbuf),(n));             \
} while (0)

void layout_swap_same_type(layout_ctx_t * restrict ctx, count_t i, count_t j)
{
    SlotA si = ctx->slots[i];
    SlotA sj = ctx->slots[j];
    assert(si.type == sj.type);

    switch (si.type) {
    case PT_GAS: {
        GasFullA tmp;
        memcpy(&tmp, &ctx->gas[si.type_idx], sizeof(tmp));
        memcpy(&ctx->gas[si.type_idx], &ctx->gas[sj.type_idx], sizeof(tmp));
        memcpy(&ctx->gas[sj.type_idx], &tmp, sizeof(tmp));
        break;
    }
    case PT_DM1: case PT_DM2: case PT_DM3: {
        DMFullA tmp;
        memcpy(&tmp, &ctx->dm[si.type_idx], sizeof(tmp));
        memcpy(&ctx->dm[si.type_idx], &ctx->dm[sj.type_idx], sizeof(tmp));
        memcpy(&ctx->dm[sj.type_idx], &tmp, sizeof(tmp));
        break;
    }
    case PT_STAR: {
        StarFullA tmp;
        memcpy(&tmp, &ctx->star[si.type_idx], sizeof(tmp));
        memcpy(&ctx->star[si.type_idx], &ctx->star[sj.type_idx], sizeof(tmp));
        memcpy(&ctx->star[sj.type_idx], &tmp, sizeof(tmp));
        break;
    }
    case PT_BH: {
        BHFullA tmp;
        memcpy(&tmp, &ctx->bh[si.type_idx], sizeof(tmp));
        memcpy(&ctx->bh[si.type_idx], &ctx->bh[sj.type_idx], sizeof(tmp));
        memcpy(&ctx->bh[sj.type_idx], &tmp, sizeof(tmp));
        break;
    }
    default:
        break;
    }
}

/* ------------------------------------------------------------------------ */
/* Accounting                                                               */
/* ------------------------------------------------------------------------ */

size_t layout_reshuffle_bytes(const layout_ctx_t *ctx)
{
    size_t payload =
        (size_t)ctx->counts.n_gas  * sizeof(GasFullA)  +
        (size_t)ctx->counts.n_dm   * sizeof(DMFullA)   +
        (size_t)ctx->counts.n_star * sizeof(StarFullA) +
        (size_t)ctx->counts.n_bh   * sizeof(BHFullA);
    size_t slot_copy = (size_t)total_parts(&ctx->counts) * sizeof(SlotA);
    return og3_perm_dense_copy_bytes(1, payload) + slot_copy;
}

int layout_reshuffle_streams(const layout_ctx_t *ctx)
{
    int n = 0;
    if (ctx->counts.n_gas)  n++;
    if (ctx->counts.n_dm)   n++;
    if (ctx->counts.n_star) n++;
    if (ctx->counts.n_bh)   n++;
    return n;
}

size_t layout_swap_bytes(const layout_ctx_t *ctx, uint8_t t)
{
    (void)ctx;
    size_t s = 0;
    switch (t) {
        case PT_GAS: s = sizeof(GasFullA); break;
        case PT_DM1: case PT_DM2: case PT_DM3: s = sizeof(DMFullA); break;
        case PT_STAR: s = sizeof(StarFullA); break;
        case PT_BH:   s = sizeof(BHFullA);   break;
    }
    /* One swap = 3 reads + 3 writes of the struct => 6 * sizeof. */
    return 6u * s;
}

int layout_swap_streams(const layout_ctx_t *ctx, uint8_t t)
{
    (void)ctx; (void)t;
    return 1;
}

/* ------------------------------------------------------------------------ */
/* Verification accessors                                                   */
/* ------------------------------------------------------------------------ */

static const GasFullA *slot_gas(const layout_ctx_t *ctx, SlotA s) { return &ctx->gas[s.type_idx]; }
static const DMFullA *slot_dm(const layout_ctx_t *ctx, SlotA s) { return &ctx->dm[s.type_idx]; }
static const StarFullA *slot_star(const layout_ctx_t *ctx, SlotA s) { return &ctx->star[s.type_idx]; }
static const BHFullA *slot_bh(const layout_ctx_t *ctx, SlotA s) { return &ctx->bh[s.type_idx]; }

pkey_t layout_get_key(const layout_ctx_t *ctx, count_t i)
{
    SlotA s = ctx->slots[i];
    switch (s.type) {
    case PT_GAS: return slot_gas(ctx, s)->key;
    case PT_DM1: case PT_DM2: case PT_DM3: return slot_dm(ctx, s)->key;
    case PT_STAR: return slot_star(ctx, s)->key;
    case PT_BH: return slot_bh(ctx, s)->key;
    default: return 0;
    }
}

void layout_set_key(layout_ctx_t *ctx, count_t i, pkey_t key)
{
    SlotA s = ctx->slots[i];
    switch (s.type) {
    case PT_GAS: ctx->gas[s.type_idx].key = key; break;
    case PT_DM1: case PT_DM2: case PT_DM3: ctx->dm[s.type_idx].key = key; break;
    case PT_STAR: ctx->star[s.type_idx].key = key; break;
    case PT_BH: ctx->bh[s.type_idx].key = key; break;
    default: break;
    }
}

void layout_get_pos(const layout_ctx_t *ctx, count_t i, pos_t out[restrict 3])
{
    SlotA s = ctx->slots[i];
    const pos_t *pos = NULL;
    switch (s.type) {
    case PT_GAS: pos = slot_gas(ctx, s)->pos; break;
    case PT_DM1: case PT_DM2: case PT_DM3: pos = slot_dm(ctx, s)->pos; break;
    case PT_STAR: pos = slot_star(ctx, s)->pos; break;
    case PT_BH: pos = slot_bh(ctx, s)->pos; break;
    default: break;
    }
    if (pos) {
        out[0] = pos[0];
        out[1] = pos[1];
        out[2] = pos[2];
    }
}

void layout_set_pos(layout_ctx_t *ctx, count_t i, const pos_t in[restrict 3])
{
    SlotA s = ctx->slots[i];
    pos_t *pos = NULL;
    switch (s.type) {
    case PT_GAS: pos = ctx->gas[s.type_idx].pos; break;
    case PT_DM1: case PT_DM2: case PT_DM3: pos = ctx->dm[s.type_idx].pos; break;
    case PT_STAR: pos = ctx->star[s.type_idx].pos; break;
    case PT_BH: pos = ctx->bh[s.type_idx].pos; break;
    default: break;
    }
    if (pos) {
        pos[0] = in[0];
        pos[1] = in[1];
        pos[2] = in[2];
    }
}

uint8_t layout_get_type(const layout_ctx_t *ctx, count_t i)
{
    return ctx->slots[i].type;
}

ppid_t layout_get_id(const layout_ctx_t *ctx, count_t i)
{
    SlotA s = ctx->slots[i];
    switch (s.type) {
    case PT_GAS: return slot_gas(ctx, s)->id;
    case PT_DM1: case PT_DM2: case PT_DM3: return slot_dm(ctx, s)->id;
    case PT_STAR: return slot_star(ctx, s)->id;
    case PT_BH: return slot_bh(ctx, s)->id;
    default: return 0;
    }
}

int layout_verify_deep(const layout_ctx_t *ctx)
{
#ifndef DEBUG
    (void)ctx;
    return 0;
#else
    count_t N = total_parts(&ctx->counts);
    for (count_t i = 0; i < N; ++i) {
        SlotA s = ctx->slots[i];
        if (s.type == PT_GAS) {
            const GasFullA *g = slot_gas(ctx, s);
            if (g->type != PT_GAS || g->swallow_id != og3_debug_marker64(g->id, 1u)) {
                fprintf(stderr, "A deep verify failed at slot %llu\n", (unsigned long long)i);
                return 1;
            }
        } else if (PT_IS_DM(s.type)) {
            const DMFullA *d = slot_dm(ctx, s);
            if (d->type != s.type) {
                fprintf(stderr, "A deep verify failed at slot %llu\n", (unsigned long long)i);
                return 1;
            }
        } else if (s.type == PT_STAR) {
            const StarFullA *st = slot_star(ctx, s);
            if (st->type != PT_STAR || st->pid_star != og3_debug_marker32(st->id, 2u)) {
                fprintf(stderr, "A deep verify failed at slot %llu\n", (unsigned long long)i);
                return 1;
            }
        } else if (s.type == PT_BH) {
            const BHFullA *bh = slot_bh(ctx, s);
            if (bh->type != PT_BH || bh->swallow_id_bh != og3_debug_marker64(bh->id, 3u)) {
                fprintf(stderr, "A deep verify failed at slot %llu\n", (unsigned long long)i);
                return 1;
            }
        }
    }
    return 0;
#endif
}
