CC ?= cc
CFLAGS ?= -std=gnu11 -O3 -march=native -Wall -Wextra -fno-strict-aliasing
LDFLAGS ?=
LDLIBS ?= -lm
LAYOUTS := A Bp Bc C
COMMON_H := bench_common.h layout_api.h fields.h hilbert.h particle_init.h boxleaves.h
RESHUFFLE_BINS := $(foreach l,$(LAYOUTS),bench_reshuffle_$(l))
PARTITION_BINS := $(foreach l,$(LAYOUTS),bench_partition_$(l))
SIZES_BINS := $(foreach l,$(LAYOUTS),bench_sizes_$(l))
ALL_BINS := $(RESHUFFLE_BINS) $(PARTITION_BINS) $(SIZES_BINS)
.PHONY: all clean info
all: $(ALL_BINS)
info:
	@echo "$(ALL_BINS)"
bench_common.o: bench_common.c $(COMMON_H)
	$(CC) $(CFLAGS) -c $< -o $@
hilbert.o: hilbert.c hilbert.h bench_common.h
	$(CC) $(CFLAGS) -c $< -o $@
particle_init.o: particle_init.c $(COMMON_H)
	$(CC) $(CFLAGS) -c $< -o $@
boxleaves.o: boxleaves.c $(COMMON_H)
	$(CC) $(CFLAGS) -c $< -o $@
layout_A.o: layout_A.c $(COMMON_H)
	$(CC) $(CFLAGS) -c $< -o $@
layout_Bp.o: layout_Bp.c $(COMMON_H)
	$(CC) $(CFLAGS) -c $< -o $@
layout_Bc.o: layout_Bc.c $(COMMON_H)
	$(CC) $(CFLAGS) -c $< -o $@
layout_C.o: layout_C.c $(COMMON_H)
	$(CC) $(CFLAGS) -c $< -o $@
bench_reshuffle_%.o: bench_reshuffle.c $(COMMON_H)
	$(CC) $(CFLAGS) -c bench_reshuffle.c -o $@
bench_partition_%.o: bench_partition.c $(COMMON_H)
	$(CC) $(CFLAGS) -c bench_partition.c -o $@
bench_sizes_%.o: bench_sizes.c $(COMMON_H)
	$(CC) $(CFLAGS) -c bench_sizes.c -o $@
bench_reshuffle_%: bench_reshuffle_%.o bench_common.o hilbert.o particle_init.o boxleaves.o layout_%.o
	$(CC) $(CFLAGS) $(LDFLAGS) $^ -o $@ $(LDLIBS)
bench_partition_%: bench_partition_%.o bench_common.o hilbert.o particle_init.o boxleaves.o layout_%.o
	$(CC) $(CFLAGS) $(LDFLAGS) $^ -o $@ $(LDLIBS)
bench_sizes_%: bench_sizes_%.o bench_common.o hilbert.o particle_init.o boxleaves.o layout_%.o
	$(CC) $(CFLAGS) $(LDFLAGS) $^ -o $@ $(LDLIBS)
clean:
	rm -f *.o $(ALL_BINS) bench_sizes
