# KH bond-direction input for honeycomb

## Goal
Define and implement input/data-layer support so every honeycomb magnetic bond can be classified as x/y/z by spatial direction (with PBC-aware supercell construction), and KH couplings can be loaded from `input.toml`, enabling downstream KH Hamiltonian terms.

## Requirements
- Target lattice for first version is honeycomb.
- POSCAR remains position-only (no bond-type labels in POSCAR itself).
- `input.toml` must add KH bond-direction configuration for x/y/z classification.
- When user selects Kitaev-Heisenberg model, `input.toml` must accept one global, flat KH coupling set `J`, `K`, `G`, `Gp` (for downstream use).
- KH coupling fields are optional and default to `0` when omitted.
- Bond direction assignment must be based on bond spatial direction, not only atom-pair identity.
- Supercell bond construction must handle periodic boundary conditions when computing bond direction.
- Internal bond object/edge record must persist direction label (`x|y|z`) for downstream usage.
- Parser must load and store one global, flat KH coupling set `J/K/G/Gp` in internal parameter structures, applying default `0` for omitted fields.
- Add validation rules for missing/invalid/ambiguous direction classification and invalid KH parameter input.
- Keep this task focused on input/data layer only (no Hamiltonian term calculations here).

## Acceptance Criteria
- [ ] Input schema supports KH bond-direction configuration and x/y/z labels.
- [ ] When KH model is selected, input schema supports one global, flat `J/K/G/Gp` set; parser loads it with default `0` for omitted fields.
- [ ] Parser loads bond-direction configuration into internal structures correctly.
- [ ] Honeycomb bonds in supercell are deterministically labeled x/y/z with PBC considered.
- [ ] Missing/invalid/ambiguous labels, classification failures, or invalid KH parameters are rejected with clear error messages.
- [ ] A minimal honeycomb example demonstrates x/y/z assignment from geometry + PBC path and includes one global KH parameter set (`J/K/G/Gp`).

## Constraints
- Spin type: classical spin (global model constraint, referenced by later tasks).
- Do not start equation research in this task.

## Task Decomposition
- Subtask A: [KH input schema for global J K G Gp](../03-31-kh-input-schema-global/prd.md)
  - Focus: flat global KH coupling fields `J/K/G/Gp`, defaults, parser mapping, parameter validation.
- Subtask B: [KH bond direction classification with PBC](../03-31-kh-bond-classify-pbc/prd.md)
  - Focus: geometry-based `x/y/z` labeling with PBC-aware bond vectors and ambiguity validation.

## Out of Scope
- K/J/Γ/Γ' energy calculations (but input parsing/storage for KH parameters is in scope).
- XSF visualization output.
