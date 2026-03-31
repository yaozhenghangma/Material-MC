# KH bond direction classification with PBC

## Goal
Implement honeycomb bond direction classification (`x|y|z`) from spatial geometry with periodic-boundary-aware supercell bond construction, and persist labels in internal bond data.

## Requirements
- Scope is honeycomb lattice first.
- POSCAR remains position-only; no manual bond labels in POSCAR.
- Bond direction assignment is computed from bond spatial direction, not atom pair identity alone.
- Bond vector construction must account for periodic boundary conditions (PBC) for supercell/neighbor bonds.
- Every relevant magnetic bond is deterministically labeled as `x`, `y`, or `z`.
- Internal bond data structure stores the resulting direction label.
- Validation reports clear errors when classification is missing/invalid/ambiguous.

## Acceptance Criteria
- [ ] Honeycomb bonds are labeled `x/y/z` deterministically in supercell context.
- [ ] PBC-crossing bonds are classified consistently with in-cell equivalent bonds.
- [ ] Direction label is persisted in internal bond records.
- [ ] Ambiguous/invalid classification triggers clear error messages.
- [ ] Minimal example demonstrates geometry+PBC classification path.

## Constraints
- Keep this task independent from KH coupling parsing logic where possible.
- Do not implement KH energy terms.

## Out of Scope
- Parsing `J/K/G/Gp` input fields.
- KH Hamiltonian term implementation.
- XSF visualization output.
