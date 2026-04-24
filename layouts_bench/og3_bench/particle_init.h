#ifndef PARTICLE_INIT_H
#define PARTICLE_INIT_H

#include "bench_common.h"
#include "layout_api.h"

#define OG3_PH_BITS 21u
#define OG3_PH_CELLS (1u << OG3_PH_BITS)

typedef struct {
    uint32_t x, y, z;
} coord3u_t;

typedef enum {
    OG3_POS_PLAIN = 0,
    OG3_POS_CLUSTERED = 1
} og3_pos_mode_t;

void og3_generate_positions(coord3u_t * restrict pos,
                            count_t N,
                            og3_pos_mode_t mode,
                            double temperature,
                            uint64_t seed);

void og3_keys_from_positions(const coord3u_t * restrict pos,
                             pkey_t * restrict keys,
                             count_t N);

void og3_set_layout_positions(layout_ctx_t * restrict ctx,
                              const coord3u_t * restrict pos,
                              count_t N);

#endif
