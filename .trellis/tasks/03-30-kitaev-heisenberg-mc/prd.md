# Kitaev-Heisenberg Monte Carlo with bond-direction input and XSF validation

## Goal
Enable Monte Carlo simulation for a Kitaev-Heisenberg family Hamiltonian that includes Kitaev, Heisenberg, Gamma, and Gamma' terms, with bond-direction-aware x/y/z link definition and visual verification via XSF output.

## Requirements
- Add input parameters to define bond directions (x/y/z) when constructing chemical bonds.
- Support couplings for K, J, Γ, and Γ' in the Hamiltonian.
- Build Hamiltonian terms according to bond direction and model conventions.
- Provide XSF output for geometry/parameter sanity check.
- In XSF output, mark x/y/z bonds with different colors.

## Proposed Task Decomposition
- Task A: Bond-direction data model and parameter input
  - Extend input schema/config parser for bond-type labeling.
  - Ensure internal bond data carries direction label (x/y/z).
  - Add validation for completeness/consistency of labels.
- Task B: Hamiltonian extension (K, J, Γ, Γ')
  - Define mathematical form used by this codebase for each interaction.
  - Map each term onto x/y/z bond-direction channels.
  - Integrate into Monte Carlo energy evaluation/update path.
- Task C: XSF visualization output
  - Add/extend XSF writer for bond visualization metadata.
  - Color-code x/y/z bonds distinctly.
  - Ensure output can be inspected by users for parameter sanity.
- Task D: Verification and examples
  - Add a minimal example input covering all four terms.
  - Add checks/tests or deterministic sanity outputs.
  - Document expected signatures in output.

## Acceptance Criteria
- [ ] User can input/select x/y/z bond directions through config.
- [ ] Monte Carlo run includes K, J, Γ, Γ' interactions correctly.
- [ ] XSF output is generated and visually distinguishes x/y/z bonds.
- [ ] A sample case is provided to verify parameter mapping.

## Confirmed Decisions
- First target lattice/topology: honeycomb.
- Spin type: classical spin.
- Γ and Γ' term formula/sign/index convention: to be determined via `mcp__grok_search` before implementation.
- XSF color mapping for bond directions: x=red, y=green, z=blue.

## Task Breakdown (Concrete)
1. `03-30-kh-bond-input` — Bond-direction input schema + parser + validation for honeycomb.
2. `03-30-kh-hamiltonian` — Classical-spin K/J/Γ/Γ' Hamiltonian integration with x/y/z bond mapping.
3. `03-30-kh-xsf-color` — XSF output with x/y/z bond color coding (RGB).
4. `03-30-kh-validation` — Minimal sample + reproducible sanity checks.

## Technical Notes
This is backend-only (C++ simulation code). The initial implementation scope is honeycomb only, with extension to other lattices deferred.
