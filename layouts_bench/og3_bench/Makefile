# =========================================================================
# Makefile for OG3 data-layout microbenchmarks.
#
# Builds 8 binaries: {reshuffle, partition} x {A, Bp, Bc, C}.
# Each binary is linked with exactly one layout_X.o, because all four
# layout TUs export the same symbol names (layout_alloc, ...).
# =========================================================================

CC      ?= cc
CFLAGS  ?= -std=gnu11 -O3 -march=native -Wall -Wextra -fno-strict-aliasing
LDFLAGS ?=
LDLIBS  ?=

LAYOUTS  := A Bp Bc C
BENCHES  := reshuffle partition

COMMON_H := bench_common.h layout_api.h fields.h hilbert.h

# Final binaries.
RESHUFFLE_BINS := $(foreach l,$(LAYOUTS),bench_reshuffle_$(l))
PARTITION_BINS := $(foreach l,$(LAYOUTS),bench_partition_$(l))
SIZES_BINS     := $(foreach l,$(LAYOUTS),bench_sizes_$(l))
ALL_BINS       := $(RESHUFFLE_BINS) $(PARTITION_BINS) $(SIZES_BINS)

.PHONY: all clean info
all: $(ALL_BINS)

info:
	@echo "Layouts  : $(LAYOUTS)"
	@echo "Benches  : $(BENCHES)"
	@echo "Binaries : $(ALL_BINS)"

# --- Common object ---
bench_common.o: bench_common.c $(COMMON_H)
	$(CC) $(CFLAGS) -c $< -o $@

hilbert.o: hilbert.c hilbert.h bench_common.h
	$(CC) $(CFLAGS) -c $< -o $@

# --- Layout objects ---
layout_A.o: layout_A.c $(COMMON_H)
	$(CC) $(CFLAGS) -c $< -o $@
layout_Bp.o: layout_Bp.c $(COMMON_H)
	$(CC) $(CFLAGS) -c $< -o $@
layout_Bc.o: layout_Bc.c $(COMMON_H)
	$(CC) $(CFLAGS) -c $< -o $@
layout_C.o:  layout_C.c  $(COMMON_H)
	$(CC) $(CFLAGS) -c $< -o $@

# --- Per-binary bench objects (one .o per binary so they don't collide) ---
bench_reshuffle_%.o: bench_reshuffle.c $(COMMON_H)
	$(CC) $(CFLAGS) -c bench_reshuffle.c -o $@
bench_partition_%.o: bench_partition.c $(COMMON_H)
	$(CC) $(CFLAGS) -c bench_partition.c -o $@
bench_sizes_%.o: bench_sizes.c $(COMMON_H)
	$(CC) $(CFLAGS) -c bench_sizes.c -o $@

# --- Link rules ---
bench_reshuffle_%: bench_reshuffle_%.o bench_common.o hilbert.o layout_%.o
	$(CC) $(CFLAGS) $(LDFLAGS) $^ -o $@ $(LDLIBS)

bench_partition_%: bench_partition_%.o bench_common.o hilbert.o layout_%.o
	$(CC) $(CFLAGS) $(LDFLAGS) $^ -o $@ $(LDLIBS)

bench_sizes_%: bench_sizes_%.o bench_common.o hilbert.o layout_%.o
	$(CC) $(CFLAGS) $(LDFLAGS) $^ -o $@ $(LDLIBS)

clean:
	rm -f *.o $(ALL_BINS) bench_sizes
