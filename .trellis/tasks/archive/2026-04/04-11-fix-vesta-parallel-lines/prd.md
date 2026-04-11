# Task: fix-vesta-parallel-lines

## Overview
Fix VESTA export behavior in `src/spin_out.cpp` so that bond segments that are geometrically parallel in the simulation data are also rendered as parallel in VESTA.

Known incorrect-output causes to address:
1. Periodic boundary wrapping can make some exported bond vectors differ by ± one supercell lattice constant (wrong image selected).
2. Some bond vectors are exported with opposite sign (direction reversed) relative to expected KH bond-direction visualization.

The fix must keep existing KH bond-color semantics unchanged.

## Requirements
- Remove coordinate-frame mixing between `STRUC` and `VECTR`/`VECTT` output in `WriteVestaKhBondColor(...)`.
- Ensure bond tail positions and bond direction vectors are written in a coordinate basis that VESTA interprets consistently with the unit cell (`CELLP`) and site coordinates.
- Correct periodic-image selection for bond vectors so nearest-neighbor KH bonds are not exported as long vectors spanning across the supercell.
- Make bond-vector orientation deterministic and consistent with the chosen KH direction-visualization convention (avoid traversal-order dependent sign flips).
- Preserve current KH color mapping in VESTA output:
  - x bond: red
  - y bond: green
  - z bond: blue
- Keep output deterministic for the same input configuration.
- Keep writer scope limited to KH VESTA output behavior (no unrelated format changes).

## Acceptance Criteria
- [ ] For a KH configuration with known parallel bonds, the generated `.vesta` file renders those bonds as parallel in VESTA.
- [ ] Code no longer mixes incompatible coordinate frames between `STRUC` and vector sections used to draw bonds.
- [ ] Bonds crossing periodic boundaries are drawn using the correct local image (no ± one-cell jump artifacts in bond vectors).
- [ ] For identical input, bond-vector orientation is stable and does not depend on neighbor traversal order.
- [ ] x/y/z bond color classes remain correctly distinguishable as red/green/blue in VESTA.
- [ ] VESTA output generation remains available from initialization and classical/PT output paths that currently invoke `WriteVestaKhBondColor(...)`.
- [ ] Output remains deterministic for identical input state.

## Technical Notes
- Current behavior to verify/fix: `STRUC` is written in fractional coordinates while `VECTR` currently uses Cartesian values derived from `BuildSiteCartesian(...)` (scaled by `magnify_factor`), which is likely causing VESTA interpretation mismatch.
- In `InitializeSupercell(...)`, neighbor pointers are materialized with modulo wrapping. Translation offsets used to reach periodic images are not retained in `Site` neighbor data, so `WriteVestaKhBondColor(...)` currently reconstructs vectors from wrapped indices only.
- `WriteVestaKhBondColor(...)` computes `delta = target_cart - source_cart` from wrapped coordinates, which can produce vectors shifted by a full cell vector when a bond crosses PBC.
- Bond de-duplication uses index-pair dedup (`visited_bonds`) rather than geometry dedup; retain deterministic representative selection.
- Current de-dup keeps the first encountered orientation for an undirected pair; this can produce sign-reversed vectors depending on iteration order.
- `WriteVestaKhBondColor(...)` only applies to `ModelType::Kitaev_Heisenberg`; reproduction cases must use KH model.
- Keep output-layer changes in `src/spin_out.cpp` and only touch upstream data generation if required for correctness.

### Brainstorm Open Decision
- PBC vector policy (confirmed): compute `d = target_frac - source_frac`, then apply periodic wrap to choose the physically shortest periodic image before converting to Cartesian.
  - In current coordinate convention this means wrapping by supercell periods (`n_x`, `n_y`, `n_z`) along each fractional axis.
  - Equivalent conceptual description: add/subtract lattice constants to obtain the shortest bond vector under PBC.
- Orientation convention for exported KH bond vectors still needs to be fixed explicitly (e.g., canonical per bond-type axis vs deterministic index-order convention).

## Out of Scope
- Hamiltonian or Monte Carlo numerical-physics changes.
- New visualization formats or UI tooling beyond existing VESTA/XSF paths.
- Broad refactoring unrelated to VESTA parallel-line rendering correctness.
