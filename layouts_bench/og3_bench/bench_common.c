/* =========================================================================
 * bench_common.c
 * ========================================================================= */

#include "bench_common.h"
#include "hilbert.h"

#include <assert.h>
#include <string.h>

void counts_from_total(count_t N_total, pcount_t * restrict out)
{
    count_t n_dm   = (N_total * 4949ULL) / 10000ULL;
    count_t n_star = (N_total *   50ULL) / 10000ULL;
    count_t n_bh   = (N_total *    1ULL) / 10000ULL;
    count_t n_gas  = N_total - n_dm - n_star - n_bh;
    out->n_gas  = n_gas;
    out->n_dm   = n_dm;
    out->n_star = n_star;
    out->n_bh   = n_bh;
}

static inline uint64_t rotl_(uint64_t x, int k)
{
    return (x << k) | (x >> (64 - k));
}

static uint64_t splitmix64(uint64_t * restrict state)
{
    uint64_t z = (*state += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

void rng_seed(rng_t * restrict r, uint64_t seed)
{
    uint64_t sm = seed;
    r->s[0] = splitmix64(&sm);
    r->s[1] = splitmix64(&sm);
    r->s[2] = splitmix64(&sm);
    r->s[3] = splitmix64(&sm);
    if ((r->s[0] | r->s[1] | r->s[2] | r->s[3]) == 0) r->s[0] = 1;
}

uint64_t rng_next(rng_t * restrict r)
{
    const uint64_t result = rotl_(r->s[1] * 5, 7) * 9;
    const uint64_t t      = r->s[1] << 17;
    r->s[2] ^= r->s[0];
    r->s[3] ^= r->s[1];
    r->s[1] ^= r->s[2];
    r->s[0] ^= r->s[3];
    r->s[2] ^= t;
    r->s[3]  = rotl_(r->s[3], 45);
    return result;
}

double rng_uniform(rng_t * restrict r)
{
    return (double)(rng_next(r) >> 11) * (1.0 / 9007199254740992.0);
}

pkey_t ph_key_from_ijk(uint32_t i, uint32_t j, uint32_t k)
{
    return hilbert_key_3d_ijk(i, j, k, 21u);
}

void gen_keys(pkey_t * restrict keys, count_t N, uint32_t n_cells, uint64_t seed)
{
    if (n_cells == 0) n_cells = 1;
    if (n_cells > (1u << 21)) n_cells = (1u << 21);

    rng_t r;
    rng_seed(&r, seed);
    for (count_t p = 0; p < N; ++p) {
        uint32_t i = (uint32_t)(rng_next(&r) % n_cells);
        uint32_t j = (uint32_t)(rng_next(&r) % n_cells);
        uint32_t k = (uint32_t)(rng_next(&r) % n_cells);
        keys[p] = ph_key_from_ijk(i, j, k);
    }
}

void assign_types(uint8_t * restrict types,
                  const pcount_t * restrict c, uint64_t seed)
{
    count_t N = total_parts(c);
    count_t dm_each = c->n_dm / 3;
    count_t dm_res  = c->n_dm - 2 * dm_each;

    count_t p = 0;
    for (count_t i = 0; i < c->n_gas;  ++i) types[p++] = PT_GAS;
    for (count_t i = 0; i < dm_each;   ++i) types[p++] = PT_DM1;
    for (count_t i = 0; i < dm_each;   ++i) types[p++] = PT_DM2;
    for (count_t i = 0; i < dm_res;    ++i) types[p++] = PT_DM3;
    for (count_t i = 0; i < c->n_star; ++i) types[p++] = PT_STAR;
    for (count_t i = 0; i < c->n_bh;   ++i) types[p++] = PT_BH;
    assert(p == N);

    if (N < 2) return;

    rng_t r;
    rng_seed(&r, seed);
    for (count_t i = N - 1; i > 0; --i) {
        count_t j = rng_next(&r) % (i + 1);
        uint8_t tmp = types[i];
        types[i] = types[j];
        types[j] = tmp;
    }
}

void radix_sort_u64_count(const pkey_t  * restrict keys_in,
                          const count_t * restrict pl_in,
                          pkey_t        * restrict keys_out,
                          count_t       * restrict pl_out,
                          count_t N)
{
    pkey_t  *k_tmp = (pkey_t  *)og3_aligned_alloc((size_t)N * sizeof(pkey_t));
    count_t *p_tmp = (count_t *)og3_aligned_alloc((size_t)N * sizeof(count_t));
    if ((N && (!k_tmp || !p_tmp))) {
        fprintf(stderr, "radix_sort: OOM\n");
        exit(1);
    }

    enum { BITS = 16, BUCKETS = 1u << BITS };
    count_t counts[BUCKETS];

    const pkey_t  *src_k = keys_in;
    const count_t *src_p = pl_in;
    pkey_t        *dst_k = NULL;
    count_t       *dst_p = NULL;

    for (int pass = 0; pass < 4; ++pass) {
        const int shift = pass * BITS;
        int to_tmp = ((pass & 1) == 0);
        dst_k = to_tmp ? k_tmp    : keys_out;
        dst_p = to_tmp ? p_tmp    : pl_out;

        memset(counts, 0, sizeof(counts));
        for (count_t i = 0; i < N; ++i) {
            uint32_t bucket = (uint32_t)((src_k[i] >> shift) & (BUCKETS - 1));
            counts[bucket]++;
        }

        count_t running = 0;
        for (uint32_t b = 0; b < BUCKETS; ++b) {
            count_t c = counts[b];
            counts[b] = running;
            running += c;
        }

        for (count_t i = 0; i < N; ++i) {
            uint32_t bucket = (uint32_t)((src_k[i] >> shift) & (BUCKETS - 1));
            count_t pos = counts[bucket]++;
            dst_k[pos] = src_k[i];
            dst_p[pos] = src_p[i];
        }

        src_k = dst_k;
        src_p = dst_p;
    }

    free(k_tmp);
    free(p_tmp);
}

void invert_permutation(const count_t * restrict perm,
                        count_t       * restrict inv,
                        count_t N)
{
    for (count_t i = 0; i < N; ++i) inv[perm[i]] = i;
}

int verify_sort(const pkey_t * restrict keys_sorted, count_t N)
{
    for (count_t i = 1; i < N; ++i) {
        if (keys_sorted[i - 1] > keys_sorted[i]) {
            fprintf(stderr, "verify_sort: out-of-order at i=%llu (%llu > %llu)\n",
                    (unsigned long long)i,
                    (unsigned long long)keys_sorted[i - 1],
                    (unsigned long long)keys_sorted[i]);
            return 1;
        }
    }
    return 0;
}

void gen_sparse_permutation(count_t * restrict perm, count_t N,
                            double fraction, uint64_t seed)
{
    for (count_t i = 0; i < N; ++i) perm[i] = i;
    if (fraction <= 0.0 || N < 2) return;

    count_t k = (count_t)(fraction * (double)N);
    if (k < 2) return;
    if (k > N) k = N;

    rng_t r;
    rng_seed(&r, seed);

    count_t *selected = (count_t *)og3_aligned_alloc((size_t)k * sizeof(count_t));
    if (!selected) { fprintf(stderr, "gen_sparse_perm: OOM\n"); exit(1); }

    for (count_t i = 0; i < k; ++i) selected[i] = i;
    for (count_t i = k; i < N; ++i) {
        count_t j = rng_next(&r) % (i + 1);
        if (j < k) selected[j] = i;
    }

    for (count_t i = k - 1; i > 0; --i) {
        count_t j = rng_next(&r) % (i + 1);
        count_t a = perm[selected[i]];
        count_t b = perm[selected[j]];
        perm[selected[i]] = b;
        perm[selected[j]] = a;
    }

    free(selected);
}

void gen_swap_pairs_same_type(swap_pair_t   * restrict pairs,
                              count_t n_pairs,
                              const uint8_t * restrict types,
                              count_t N, uint64_t seed)
{
    count_t count[NTYPES] = {0};
    for (count_t i = 0; i < N; ++i) count[types[i]]++;

    count_t *idx_store = (count_t *)og3_aligned_alloc((size_t)N * sizeof(count_t));
    if (N && !idx_store) { fprintf(stderr, "gen_swap_pairs: OOM\n"); exit(1); }

    count_t *idx[NTYPES];
    count_t off = 0;
    for (int t = 0; t < NTYPES; ++t) {
        idx[t] = idx_store + off;
        off += count[t];
    }

    count_t fill[NTYPES] = {0};
    for (count_t i = 0; i < N; ++i) {
        uint8_t t = types[i];
        idx[t][fill[t]++] = i;
    }

    count_t pool_total = 0;
    count_t cum[NTYPES];
    for (int t = 0; t < NTYPES; ++t) {
        if (count[t] >= 2) pool_total += count[t];
        cum[t] = pool_total;
    }
    if (pool_total == 0) {
        fprintf(stderr, "gen_swap_pairs: no type has >=2 members\n");
        exit(1);
    }

    rng_t r;
    rng_seed(&r, seed);
    for (count_t p = 0; p < n_pairs; ++p) {
        count_t pick = rng_next(&r) % pool_total;
        int t = 0;
        while (t < NTYPES && cum[t] <= pick) ++t;

        count_t i = idx[t][rng_next(&r) % count[t]];
        count_t j;
        do {
            j = idx[t][rng_next(&r) % count[t]];
        } while (j == i);
        pairs[p].i = i;
        pairs[p].j = j;
    }

    free(idx_store);
}

count_t count_permutation_changed(const count_t * restrict perm, count_t N)
{
    count_t changed = 0;
    for (count_t i = 0; i < N; ++i) changed += (perm[i] != i);
    return changed;
}

void og3_apply_perm_bytes(void * restrict array,
                          const count_t * restrict perm,
                          count_t N,
                          size_t elem_size,
                          void * restrict tmp_bytes,
                          count_t * restrict changed_idx)
{
    char *src = (char *)array;
    char *tmp = (char *)tmp_bytes;

    count_t changed = 0;
    for (count_t i = 0; i < N; ++i) {
        if (perm[i] != i) changed_idx[changed++] = i;
    }
    if (changed == 0) return;

    /* Dense path once most of the array moves. This avoids an extra pass over
     * changed_idx and gives the compiler a simple streaming copy loop. */
    if (changed * 4 >= N) {
        for (count_t i = 0; i < N; ++i) {
            memcpy(tmp + (size_t)i * elem_size,
                   src + (size_t)perm[i] * elem_size,
                   elem_size);
        }
        memcpy(src, tmp, (size_t)N * elem_size);
        return;
    }

    for (count_t m = 0; m < changed; ++m) {
        count_t dst = changed_idx[m];
        memcpy(tmp + (size_t)m * elem_size,
               src + (size_t)perm[dst] * elem_size,
               elem_size);
    }
    for (count_t m = 0; m < changed; ++m) {
        count_t dst = changed_idx[m];
        memcpy(src + (size_t)dst * elem_size,
               tmp + (size_t)m * elem_size,
               elem_size);
    }
}

size_t og3_perm_dense_copy_bytes(count_t N, size_t elem_size)
{
    return 4ull * (size_t)N * elem_size;
}

size_t og3_perm_sparse_copy_bytes(count_t changed, size_t elem_size)
{
    return 2ull * (size_t)changed * elem_size;
}

#ifdef DEBUG
uint64_t og3_debug_marker64(ppid_t id, uint64_t tag)
{
    uint64_t x = (uint64_t)id ^ (0x9E3779B97F4A7C15ULL + (tag << 32));
    x ^= x >> 33;
    x *= 0xff51afd7ed558ccdULL;
    x ^= x >> 33;
    return x;
}

uint32_t og3_debug_marker32(ppid_t id, uint32_t tag)
{
    uint64_t x = og3_debug_marker64(id, tag);
    return (uint32_t)(x ^ (x >> 32));
}
#endif
