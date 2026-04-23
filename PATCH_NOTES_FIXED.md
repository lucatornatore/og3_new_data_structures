# Patch notes for the fixed prototype

This patch set focuses on correctness, resident-layer gating, and arena ergonomics.

## What changed

- `opg_common/permutation/permutation.hpp`
  - reshuffle now gathers through byte scratch (`void*` + `memcpy`), so caller scratch alignment is no longer a correctness hazard;
  - element swaps are byte-wise too, preserving padding deterministically;
  - `ArrayRegistry` gained `swap_all(a, b)` so mutation paths can move every registered array automatically.

- `opg_common/memory/memory_arena.hpp`
  - the arena now supports two disciplines in one buffer:
    - persistent bottom-up allocations that can be individually returned and later reused;
    - top-of-buffer stack allocations for temporary work buffers via `allocate_temp*()` + `ArenaCheckpoint`.

- `opg_common/types/physics_config.hpp`
  - `stellar_evolution` and `star_formation` are treated as a compile-time invariant;
  - element-abundance storage is enabled only when `n_metal_species > 0`;
  - cosmic-ray storage is enabled only when the feature is on and at least one bin count is non-zero.

- `opg_layout_b`, `opg_layout_c`, `opg_layout_b_coarse`
  - gas-side resident storage is now absent from the API and from allocation when hydro is disabled;
  - `NMET=0` now gives energy-only stellar evolution without element arrays in `StarCore<0>` / `GasMetal<0>`;
  - coarse-B regained deterministic global reshuffle by using the shared PH-key sorter;
  - coarse-B mutation now recomputes all leaf offsets after mutation and repairs common-array order automatically;
  - C mutation now uses registry-driven automatic swaps for both common and gas registries.

- examples
  - B and coarse-B examples now allocate reshuffle scratch from the arena temp stack.

## Regression checks run

- full project build with CMake;
- shipped examples:
  - `example_layout_b`
  - `example_layout_c`
  - `example_layout_b_coarse`
- negative compile tests showing `gas_core()` / `gas_all()` are absent for `GravityOnly`;
- `NMET=0` size check showing `StarCore<0>` is smaller than `StarCore<11>` in all three layouts;
- coarse-B mutation check verifying repaired leaf ordering and offsets;
- arena check verifying free-list reuse plus temp-stack restore.
