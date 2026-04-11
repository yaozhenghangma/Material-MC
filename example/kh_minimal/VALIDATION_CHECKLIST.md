# KH minimal validation checklist

Use this checklist to sanity-check KH parsing, KH Hamiltonian runtime dispatch,
and x/y/z bond-color output with a lightweight reproducible input.

## Sample artifacts

- Structure: `example/kh_minimal/POSCAR`
- Input: `example/kh_minimal/input.toml`

The input intentionally covers all KH couplings (`J`, `K`, `G`, `Gp`) and a
full bond-direction map (`type1/type2/type3 -> x/y/z`).

## Run command

From project root:

```bash
mpirun -np 1 ./build/MMC -c example/kh_minimal/POSCAR -i example/kh_minimal/input.toml -o kh_validation_output.txt -s kh_validation_spin
```

## Checkpoint A — bond-direction parsing sanity

Check `log.txt` contains:

- `Model Type: Kitaev-Heisenberg model.`
- `KH global couplings (J, K, G, Gp):     -1.00000      2.00000      0.50000     -0.25000`
- `KH bond type mapping (type1, type2, type3): x y z`

## Checkpoint B — Hamiltonian path invocation sanity

Check `log.txt` contains:

- `Hamiltonian: Kitaev-Heisenberg model (bond-dependent J/K/G/Gp terms).`

This confirms KH Hamiltonian dispatch is active.

## Checkpoint C — x/y/z color output sanity

Check output file exists:

- `structure_initialized_kh_bond_color.vesta`

In its `VECTT` section, verify color coding:

- x-direction bonds: `255 0 0` (red)
- y-direction bonds: `0 255 0` (green)
- z-direction bonds: `0 0 255` (blue)

## Reproducibility checkpoint

Run the same command twice and compare checksums:

```bash
shasum log.txt structure_initialized_kh_bond_color.vesta
```

For identical inputs/environment, checksums should be identical.
