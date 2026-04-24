/* =========================================================================
 * layout_api.h
 *
 * Uniform API exported by layout_A.c, layout_Bp.c, layout_Bc.c, and
 * layout_C.c. The permutation convention is gather-style:
 *   after reshuffle, new slot i contains old slot perm[i].
 * ========================================================================= */

#ifndef LAYOUT_API_H
#define LAYOUT_API_H

#include "bench_common.h"

#include <stddef.h>
#include <stdint.h>

typedef struct layout_ctx layout_ctx_t;

const char *layout_name(void);
const char *layout_description(void);

layout_ctx_t *layout_alloc(const pcount_t * restrict c);
void layout_free(layout_ctx_t * restrict ctx);

void layout_fill(layout_ctx_t * restrict ctx,
                 const pkey_t  * restrict keys,
                 const uint8_t * restrict types,
                 const ppid_t  * restrict ids);

void layout_reshuffle_full(layout_ctx_t * restrict ctx,
                           const count_t * restrict perm);

void layout_swap_same_type(layout_ctx_t * restrict ctx, count_t i, count_t j);

/* Estimated actual copy traffic for the implementation's dense reshuffle path,
 * including gather read, scratch write, scratch read, and final write for all
 * permuted payload streams. Index-array traffic is reported by the benchmark. */
size_t layout_reshuffle_bytes(const layout_ctx_t * restrict ctx);
int    layout_reshuffle_streams(const layout_ctx_t * restrict ctx);

size_t layout_swap_bytes(const layout_ctx_t * restrict ctx, uint8_t t);
int    layout_swap_streams(const layout_ctx_t * restrict ctx, uint8_t t);

pkey_t  layout_get_key (const layout_ctx_t * restrict ctx, count_t i);
void    layout_set_key (layout_ctx_t * restrict ctx, count_t i, pkey_t key);

void    layout_get_pos (const layout_ctx_t * restrict ctx, count_t i,
                        pos_t out[restrict 3]);
void    layout_set_pos (layout_ctx_t * restrict ctx, count_t i,
                        const pos_t in[restrict 3]);

uint8_t layout_get_type(const layout_ctx_t * restrict ctx, count_t i);
ppid_t  layout_get_id  (const layout_ctx_t * restrict ctx, count_t i);

/* Expensive semantic verification. In non-DEBUG builds this is a zero-cost
 * stub returning success; in DEBUG builds it checks cross-index/positional
 * markers in the type-specific payloads. */
int layout_verify_deep(const layout_ctx_t * restrict ctx);

#endif /* LAYOUT_API_H */
