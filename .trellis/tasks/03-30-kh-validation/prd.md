# KH validation sample and checks

## Goal
Add validation artifacts that let users verify KH parameter setup (bond directions + K/J/Γ/Γ') and corresponding outputs.

## Requirements
- Provide at least one minimal honeycomb sample configuration covering K/J/Γ/Γ'.
- Include expected checkpoints for:
  - bond-direction parsing sanity
  - Hamiltonian path invocation sanity
  - XSF color output sanity (red/green/blue for x/y/z)
- Keep checks lightweight and reproducible.

## Acceptance Criteria
- [ ] Sample input(s) exist and run through the target flow.
- [ ] Validation checklist is documented in task artifact.
- [ ] Re-running produces consistent outputs for identical input.

## Dependencies
- Depends on: `03-30-kh-bond-input`, `03-30-kh-hamiltonian`, `03-30-kh-xsf-color`.

## Out of Scope
- Large-scale benchmarking.
- Full physics verification against external published phase diagrams.
