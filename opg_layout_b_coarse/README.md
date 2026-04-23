# opg_layout_b_coarse - Layout B (coarse variant)

A coarse-grained counterpart to `opg_layout_b` (B'). It keeps the same
explicit cross-index scheme via `PLinkageB`, but groups resident data into a
small number of wider arrays.

## Why this variant exists

The testbed needs to disentangle two axes independently:
- connection scheme: cross-index (B/B') versus positional mapping (C)
- granularity: coarse (B) versus fine (B'/C)

Layout B-coarse isolates the granularity question while keeping the simpler
cross-index mutation story.

## Resident layout

Common side:
- `PCoreB`      : 64 B
- `PDynB<Cfg>`  : 64 B
- `PAuxB<Cfg>`  : 64 B
- `PLinkageB`   : 8 B, local-only, never sent on the wire

Type-specific side:
- `GasAllB<Cfg>`
- `StarAllB<Cfg>`
- `BHAllB<Cfg>`

PH-key reshuffle touches only the four common arrays. The typed arrays are not
permuted; linkage rides with the common slot.

## Tight MPI bridge

The current bridge sends one contiguous `MPI_BYTE` message per neighbour.

The wire format keeps a single dense common block. Each common record is the
byte concatenation of `(PCoreB, PDynB<Cfg>, PAuxB<Cfg>)`, repeated once per
outgoing particle.

```text
[ PackedExchangeHeader ]
[ common block = (PCoreB, PDynB<Cfg>, PAuxB<Cfg>) repeated for each outgoing particle ]
[ GasAllB<Cfg>[n_gas] ]
[ StarAllB<Cfg>[n_star] ]
[ BHAllB<Cfg>[n_bh] ]
```

This removes the old fixed-record waste where DM particles still carried a
max-sized typed payload on the wire.

### Sender side

1. Count outgoing particles by type.
2. Compute exact bytes with `exchange_bytes<Cfg>(counts)`.
3. Allocate one temp buffer from `MemoryArena`.
4. Call `pack_range_exchange(...)` or `pack_scatter_exchange(...)`.
5. Send the resulting `(void*, bytes)` as `MPI_BYTE`.

### Receiver side

1. Receive one byte buffer.
2. Load and validate `PackedExchangeHeader`.
3. Call `unpack_exchange(...)`.
4. The receiver reconstructs `PLinkageB` locally from the type tags carried in
   `PCoreB`.

## Example result

For the synthetic StandardSPH smoke test with 37 particles
(10 gas, 20 DM, 5 star, 2 BH):

- old fixed-record bridge: 18,944 B
- tight bridge: 11,024 B
- reduction: 41.8%

## Files

```text
opg_layout_b_coarse/include/opg_layout_b_coarse/
|- particles.hpp
|- particle_container.hpp
|- reshuffle.hpp
|- pack_unpack.hpp
|- mutation.hpp
`- compact.hpp
```

Utilities:
- `example_layout_b_coarse` demonstrates reshuffle plus tight pack/unpack
- `print_exchange_size` reports exact wire sizes for a chosen particle mix
