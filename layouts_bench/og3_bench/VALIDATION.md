# Validation - v5 partition rewrite

Performed in the container on the source tree in `layouts_bench/og3_bench`.

## Manual build used for smoke tests

The full Makefile remains present. For fast smoke validation, the partition
objects were compiled manually at `-O0`:

```sh
cc -std=gnu11 -O0 -g -Wall -Wextra -fno-strict-aliasing -c bench_common.c
cc -std=gnu11 -O0 -g -Wall -Wextra -fno-strict-aliasing -c hilbert.c
cc -std=gnu11 -O0 -g -Wall -Wextra -fno-strict-aliasing -c particle_init.c
cc -std=gnu11 -O0 -g -Wall -Wextra -fno-strict-aliasing -c boxleaves.c
cc -std=gnu11 -O0 -g -Wall -Wextra -fno-strict-aliasing -c bench_partition.c
cc -std=gnu11 -O0 -g -Wall -Wextra -fno-strict-aliasing -c layout_A.c
cc -std=gnu11 -O0 -g -Wall -Wextra -fno-strict-aliasing -c layout_Bp.c
cc -std=gnu11 -O0 -g -Wall -Wextra -fno-strict-aliasing -c layout_Bc.c
cc -std=gnu11 -O0 -g -Wall -Wextra -fno-strict-aliasing -c layout_C.c
```

## Smoke tests passed

Plain mode:

```sh
./bench_partition_A  1000 1 64 32 0.5 2097152 0 0.25
./bench_partition_Bp 1000 1 64 32 0.5 2097152 0 0.25
./bench_partition_Bc 1000 1 64 32 0.5 2097152 0 0.25
./bench_partition_C  1000 1 64 32 0.5 2097152 0 0.25
```

Small-N clustered mode, where `N < target`:

```sh
./bench_partition_A  20 1 2 64 0.5 2097152 1 0.05
./bench_partition_Bp 20 1 2 64 0.5 2097152 1 0.05
./bench_partition_Bc 20 1 2 64 0.5 2097152 1 0.05
./bench_partition_C  20 1 2 64 0.5 2097152 1 0.05
```

Clustered 32^3-sized test:

```sh
./bench_partition_Bc 32768 1 32 1024 0.5 2097152 1 0.05
```

Max-side stress test:

```sh
./bench_partition_Bc 1000 1 64 32 0.5 1024 0 0.25
```

All completed with `leaf partitioning + gas subranges: OK`.

## Not completed in this pass

A full `make all` optimized build and full DEBUG/ASan matrix were not completed
in this run. The changed source files compile in the manual smoke build above,
and the new partition path passed runtime semantic checks for all four layouts.
