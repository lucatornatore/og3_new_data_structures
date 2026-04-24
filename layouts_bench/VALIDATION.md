# Validation performed for v2 edit

Commands run in this environment:

```bash
make clean && make all
```

Result: all 12 binaries built with the default optimized flags.

```bash
make clean && \
CFLAGS='-std=gnu11 -O0 -g -Wall -Wextra -DDEBUG -fno-strict-aliasing' \
make all
```

Result: all 12 binaries built with DEBUG marker instrumentation enabled.

```bash
make clean && \
CFLAGS='-std=gnu11 -O2 -Wall -Wextra -Wpedantic -fno-strict-aliasing' \
make all
```

Result: all 12 binaries built with `-Wpedantic`.

```bash
make clean && \
CFLAGS='-std=gnu11 -O1 -g -Wall -Wextra -DDEBUG \
        -fsanitize=address,undefined -fno-omit-frame-pointer \
        -fno-strict-aliasing' \
LDFLAGS='-fsanitize=address,undefined' \
make all
```

Then, for every layout `A Bp Bc C`:

```bash
./bench_reshuffle_$x 2000 1 1.0
./bench_reshuffle_$x 2000 1 0.1
./bench_partition_$x 2000 0.02 1
```

Result: AddressSanitizer/UBSan smoke tests passed for dense reshuffle, sparse
reshuffle, and same-type swap primitive on all four layouts.

Additional optimized smoke tests:

```bash
./bench_reshuffle_$x 100000 3 1.0
./bench_reshuffle_$x 100000 3 0.1
./bench_partition_$x 100000 0.01 3
```

Result: semantic verification passed for all four layouts. Timings in this
sandbox are coarse because `CLOCK_PROCESS_CPUTIME_ID` appears quantized at the
small test sizes; use larger `N` and/or `n_iter` for benchmark numbers.


# Validation performed for key/position API edit

This edit added the following public API functions to `layout_api.h` and all
four layout implementations:

```c
void layout_set_key(layout_ctx_t *ctx, count_t i, pkey_t key);
void layout_get_pos(const layout_ctx_t *ctx, count_t i, pos_t out[3]);
void layout_set_pos(layout_ctx_t *ctx, count_t i, const pos_t in[3]);
```

`bench_reshuffle.c` now performs a round-trip self-test immediately after
`layout_fill()` and before the first reshuffle: it temporarily changes and
restores the key and position for a few global slots using only the public
layout API.

Commands run in this environment after the edit:

```bash
make clean && make all
```

Result: all 12 binaries built with the default optimized flags.

```bash
make clean && CFLAGS='-std=gnu11 -O0 -g -Wall -Wextra -DDEBUG -fno-strict-aliasing' make all
```

Result: all 12 binaries built with DEBUG marker instrumentation enabled.

Then, for every layout `A Bp Bc C`:

```bash
./bench_reshuffle_$x 2000 1 1.0
./bench_reshuffle_$x 2000 1 0.1
./bench_partition_$x 2000 0.02 1
```

Result: semantic checks, deep DEBUG checks, and the new key/position API
self-test passed for dense reshuffle, sparse reshuffle, and the same-type swap
primitive on all four layouts.

```bash
make clean && CFLAGS='-std=gnu11 -O2 -Wall -Wextra -Wpedantic -fno-strict-aliasing' make all
```

Result: all 12 binaries built with `-Wpedantic`.

```bash
make clean && CFLAGS='-std=gnu11 -O1 -g -Wall -Wextra -DDEBUG         -fsanitize=address,undefined -fno-omit-frame-pointer         -fno-strict-aliasing' LDFLAGS='-fsanitize=address,undefined' make all
```

Then, for every layout `A Bp Bc C`:

```bash
./bench_reshuffle_$x 2000 1 1.0
./bench_reshuffle_$x 2000 1 0.1
./bench_partition_$x 2000 0.02 1
```

Result: AddressSanitizer/UBSan smoke tests passed for all four layouts.
