# KH Hamiltonian terms K J Gamma Gamma-prime

## Goal
Construct the classical-spin honeycomb Kitaev-Heisenberg-family Hamiltonian including K, J, Γ, and Γ' terms, using bond-direction labels x/y/z.

## Requirements
- Spin model is classical spin.
- Terms to include: Kitaev (K), Heisenberg (J), Gamma (Γ), Gamma-prime (Γ').
- Bond-dependent mapping must use x/y/z bond labels from input task.
- Formula/sign/index convention for Γ and Γ' must be determined via `mcp__grok_search` and then fixed in implementation notes before coding.
- Integrate with Monte Carlo energy path without changing unrelated model behavior.

## Acceptance Criteria
- [ ] Hamiltonian evaluation includes K/J/Γ/Γ' for honeycomb.
- [ ] Bond-direction dependence is correctly applied for each bond type.
- [ ] Γ/Γ' convention source and chosen formula are recorded in task notes/implementation artifact.
- [ ] Existing non-KH workflows remain functional.

## Dependencies
- Depends on: `03-30-kh-bond-input` (bond-direction data available).

## Out of Scope
- XSF color rendering details.
- Broader non-honeycomb generalization.

## Convention Source/Decision
- Adopted KH nearest-neighbor bond Hamiltonian convention:
  - H_ij = J Si·Sj + K Si^γ Sj^γ + G (Si^α Sj^β + Si^β Sj^α)
    + Gp (Si^α Sj^γ + Si^γ Sj^α + Si^β Sj^γ + Si^γ Sj^β)
  - Cyclic mapping: γ=x => (α,β)=(y,z), γ=y => (α,β)=(z,x), γ=z => (α,β)=(x,y)
- References consulted and aligned for Γ/Γ' index/sign convention:
  - https://arxiv.org/abs/1310.7940
  - https://link.aps.org/doi/10.1103/PhysRevLett.112.077204
  - https://jeffrau.ca/assets/papers/kitaev-2014.pdf
