# KH input schema for global J K G Gp

## Goal
Define input.toml schema and parser behavior for Kitaev-Heisenberg global couplings (`J`, `K`, `G`, `Gp`) so downstream KH Hamiltonian can consume validated parameters.

## Requirements
- Scope is input/data layer only, no energy evaluation.
- When KH model is selected, `input.toml` accepts flat global fields: `J`, `K`, `G`, `Gp`.
- Missing coupling fields default to `0`.
- Parser maps these values into internal KH parameter storage.
- Validation reports clear errors for invalid numeric type/format.

## Acceptance Criteria
- [ ] Input schema clearly defines flat global fields `J/K/G/Gp`.
- [ ] Parser loads these four parameters correctly.
- [ ] Omitted fields are assigned default value `0`.
- [ ] Invalid parameter values trigger clear error messages.
- [ ] Minimal KH input example is provided.

## Constraints
- Keep compatibility with existing non-KH model inputs.
- Do not add per-bond or per-direction coupling overrides in this task.

## Out of Scope
- Bond x/y/z classification logic.
- KH Hamiltonian term implementation.
- XSF visualization output.
