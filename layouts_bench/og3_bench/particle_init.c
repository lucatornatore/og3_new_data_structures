#include "particle_init.h"
#include "layout_api.h"

#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

static uint32_t clamp_coord_i64(int64_t v)
{
    if (v < 0) return 0;
    uint32_t hi = OG3_PH_CELLS - 1u;
    if ((uint64_t)v > (uint64_t)hi) return hi;
    return (uint32_t)v;
}

static double gaussian01(rng_t * restrict r)
{
    double u1 = rng_uniform(r);
    double u2 = rng_uniform(r);
    if (u1 < 1e-300) u1 = 1e-300;
    return sqrt(-2.0 * log(u1)) * cos(6.2831853071795864769 * u2);
}

static uint32_t uniform_coord(rng_t * restrict r)
{
    return (uint32_t)(rng_next(r) & (OG3_PH_CELLS - 1u));
}

void og3_generate_positions(coord3u_t * restrict pos,
                            count_t N,
                            og3_pos_mode_t mode,
                            double temperature,
                            uint64_t seed)
{
    /* The uploaded create_boxleaves prototype generated either plain random
     * points or clustered points near randomly chosen centres. Its OpenMP loop
     * modified the loop index inside an omp-for cluster run, which can race
     * with other iterations. Here the same two modes are implemented with
     * independent per-particle RNG streams, so the helper is safe with or
     * without OpenMP. Temperature is a coldness parameter in [0,1]: 0 is very
     * cold/tightly clustered, 1 is hot/broad. */
    if (temperature < 0.0) temperature = 0.0;
    if (temperature > 1.0) temperature = 1.0;

    enum { NCL = 16 };
    coord3u_t centres[NCL];
    rng_t cr;
    rng_seed(&cr, seed ^ UINT64_C(0x9e3779b97f4a7c15));
    for (int c = 0; c < NCL; ++c) {
        /* Keep centres away from hard boundaries, as in the uploaded generator,
         * so clusters do not all get clipped in cold runs. */
        uint32_t pad = OG3_PH_CELLS / 200u;
        uint32_t span = OG3_PH_CELLS - 2u * pad;
        centres[c].x = pad + (uint32_t)(rng_next(&cr) % span);
        centres[c].y = pad + (uint32_t)(rng_next(&cr) % span);
        centres[c].z = pad + (uint32_t)(rng_next(&cr) % span);
    }

    const double min_sigma = 2.0;
    const double max_sigma = (double)OG3_PH_CELLS / 30.0;
    double sigma = min_sigma + temperature * max_sigma;

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (count_t i = 0; i < N; ++i) {
        rng_t r;
        rng_seed(&r, seed + UINT64_C(0xd1b54a32d192ed03) * (i + 1u));

        if (mode == OG3_POS_PLAIN) {
            pos[i].x = uniform_coord(&r);
            pos[i].y = uniform_coord(&r);
            pos[i].z = uniform_coord(&r);
        } else {
            /* 95% clustered and 5% background mirrors the uploaded generator's
             * late uniform tail, but without overlapping writes. */
            int use_cluster = (rng_uniform(&r) < 0.95);
            if (!use_cluster) {
                pos[i].x = uniform_coord(&r);
                pos[i].y = uniform_coord(&r);
                pos[i].z = uniform_coord(&r);
            } else {
                int c = (int)(rng_next(&r) % NCL);
                int64_t dx = (int64_t)llround(gaussian01(&r) * sigma);
                int64_t dy = (int64_t)llround(gaussian01(&r) * sigma);
                int64_t dz = (int64_t)llround(gaussian01(&r) * sigma);
                pos[i].x = clamp_coord_i64((int64_t)centres[c].x + dx);
                pos[i].y = clamp_coord_i64((int64_t)centres[c].y + dy);
                pos[i].z = clamp_coord_i64((int64_t)centres[c].z + dz);
            }
        }
    }
}

void og3_keys_from_positions(const coord3u_t * restrict pos,
                             pkey_t * restrict keys,
                             count_t N)
{
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (count_t i = 0; i < N; ++i)
        keys[i] = ph_key_from_ijk(pos[i].x, pos[i].y, pos[i].z);
}

void og3_set_layout_positions(layout_ctx_t * restrict ctx,
                              const coord3u_t * restrict pos,
                              count_t N)
{
    for (count_t i = 0; i < N; ++i) {
        pos_t p[3] = {(pos_t)pos[i].x, (pos_t)pos[i].y, (pos_t)pos[i].z};
        layout_set_pos(ctx, i, p);
    }
}
