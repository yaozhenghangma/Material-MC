# KH XSF bond coloring output

## Goal
Provide XSF-based structure/bond output to visually verify honeycomb x/y/z bond assignments with distinct colors.

## Requirements
- Generate XSF output suitable for user-side visual inspection.
- Color mapping is fixed as:
  - x bond: red
  - y bond: green
  - z bond: blue
- Output must clearly indicate bond segments and associated color class.
- Keep writer behavior deterministic for same input.

## Acceptance Criteria
- [ ] XSF file is generated for KH honeycomb configurations.
- [ ] x/y/z bonds are distinguishable as red/green/blue in intended viewer workflow.
- [ ] Output can be used to detect wrong bond-direction setup.

## Dependencies
- Depends on: `03-30-kh-bond-input` (bond labels available).

## Out of Scope
- Hamiltonian numerical correctness (handled by Hamiltonian task).
- New visualization formats beyond XSF.
