# KH bond-direction input for honeycomb

## Goal
Define and implement parameter input/data representation for honeycomb bond directions (x/y/z) so downstream Hamiltonian terms can be applied correctly.

## Requirements
- Target lattice for first version is honeycomb.
- Input format must explicitly carry bond direction labels x/y/z.
- Internal bond object/edge record must preserve direction label.
- Add validation rules for missing/invalid direction labels.
- Keep this task focused on input/data layer only (no Hamiltonian term calculations here).

## Acceptance Criteria
- [ ] Input schema supports x/y/z labels for each relevant bond.
- [ ] Parser loads labels correctly into internal structures.
- [ ] Invalid labels are rejected with clear error message.
- [ ] A minimal honeycomb input example with x/y/z assignment is provided.

## Constraints
- Spin type: classical spin (global model constraint, referenced by later tasks).
- Do not start equation research in this task.

## Out of Scope
- K/J/Γ/Γ' energy calculations.
- XSF visualization output.
