# KH XSF bond coloring output

## Goal
Provide visualization output to verify honeycomb x/y/z bond assignments with distinct colors, using VESTA as the primary colored-bond format while keeping XSF output available.

## Requirements
- Generate VESTA output suitable for user-side visual inspection of bond classes.
- Keep existing XSF output behavior available (do not remove existing writer output path).
- Color mapping in VESTA is fixed as:
  - x bond: red
  - y bond: green
  - z bond: blue
- Output must clearly indicate bond segments and associated color class.
- Keep writer behavior deterministic for same input.

## Acceptance Criteria
- [ ] VESTA file is generated for KH honeycomb configurations.
- [ ] x/y/z bonds are distinguishable as red/green/blue in VESTA.
- [ ] Existing XSF output remains available.
- [ ] Output can be used to detect wrong bond-direction setup.

## Dependencies
- Depends on: `03-30-kh-bond-input` (bond labels available).

## Out of Scope
- Hamiltonian numerical correctness (handled by Hamiltonian task).
- New visualization formats beyond XSF.
