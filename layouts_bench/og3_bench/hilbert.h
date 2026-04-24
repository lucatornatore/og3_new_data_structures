#ifndef OG3_HILBERT_H
#define OG3_HILBERT_H

#include "bench_common.h"

/* Martin Reinecke/Gadget 3D Peano-Hilbert key. x, y, z must be in
 * [0, 2^bits). bits must be in [1, 21] for a 63-bit key. */
pkey_t hilbert_key_3d_ijk(uint32_t x, uint32_t y, uint32_t z, unsigned bits);

#endif /* OG3_HILBERT_H */
