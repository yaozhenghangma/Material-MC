Architecture Overview
=====================

This project implements a Monte Carlo simulation framework for the classical Heisenberg (and optional Ising) spin models on periodic supercells. The code is organized around a small set of core data structures and a clear simulation pipeline, with MPI used to distribute temperature points and (optionally) perform replica exchange.

High-level flow
---------------

1. **Initialization and input parsing**
   - The entry point `src/MMC.cpp` reads command-line options, a TOML input file, and a POSCAR-like structure file.
   - Input is split by responsibility:
     - `input.toml` defines Monte Carlo parameters, model/hamiltonian options, lattice size, neighbor definitions, and initialization settings.
     - `POSCAR` provides the lattice vectors and atomic positions for the base cell.

2. **Data distribution (MPI)**
   - The root process parses input and then broadcasts the `Supercell` and `MonteCarlo` configuration to all ranks.

3. **Supercell construction**
   - The lattice is enlarged according to `n_x, n_y, n_z`.
   - Spins are initialized using either explicit directions or rotations derived from `Initialization` angles.
   - Neighbor lists are built for each site by matching element types and distances with tolerance.

4. **Simulation loop**
   - The framework selects one of two methods:
     - **Classical Monte Carlo**: each rank handles a subset of temperature points.
     - **Parallel Tempering**: ranks perform replica exchange at fixed intervals.
   - Each Monte Carlo step updates a randomly chosen spin using the Metropolis criterion.

5. **Output**
   - Thermodynamic observables (energy, heat capacity, magnetization, susceptibility) are collected and written to an output file.
   - Spin configurations can be exported to XSF for visualization.

Core data structures
--------------------

- **BaseSite** (`src/MC_structure.h`)
  - Stores per-species parameters read from input: element names, exchange parameters, anisotropy, magnetic field, and neighbor definitions.

- **Site** (`src/MC_structure.h`)
  - Represents a single spin site in the supercell.
  - Holds pointers back to BaseSite parameters and a list of neighbor site pointers.

- **Lattice** (`src/MC_structure.h`)
  - Lattice vectors, scaling, supercell size, Hamiltonian type, model type, and optional external field.

- **MonteCarlo** (`src/MC_structure.h`)
  - Temperature schedule, number of steps, flips per step, and method selection.

- **Supercell** (`src/MC_structure.h`)
  - Owns the 4D container of `Site` objects.
  - Stores function pointers for the active Hamiltonian and update rule.

Hamiltonian and updates
-----------------------

- Hamiltonians are defined in `src/Hamiltonian.cpp` and selected at runtime based on input.
- The update rule is selected by model type:
  - Heisenberg: random 3D spin proposals.
  - Ising: spin flips.

Parallelization model
---------------------

- **Classical Monte Carlo**: temperatures are partitioned across ranks; each rank relaxes and samples at its assigned temperatures, then results are gathered on rank 0.
- **Parallel Tempering**: replica exchange steps occur at a configurable interval using Metropolis acceptance, with temperatures swapped between neighboring ranks.

Configuration boundaries
------------------------

- `input.toml` defines physics and simulation parameters.
- `POSCAR` defines the base lattice and atom positions.
- The code uses a tolerance-based distance check to match neighbor shells.

Key entry points and files
--------------------------

- Entry point: `src/MMC.cpp`
- Input parsing: `src/configure_in.cpp`, `src/structure_in.cpp`
- Supercell construction: `src/initialization.cpp`
- Hamiltonians: `src/Hamiltonian.cpp`
- Monte Carlo methods: `src/methods/classical.cpp`, `src/methods/parallel_tempering.cpp`
- Output: `src/result_out.cpp`, `src/spin_out.cpp`

Design notes for agents
-----------------------

- Runtime behavior is driven entirely by the parsed `Supercell` and `MonteCarlo` structures.
- The selected Hamiltonian and update strategy are injected via function pointers in `Supercell`.
- Neighbor topology is derived once during initialization and reused throughout the simulation.
- MPI is used only for distributing temperature points and replica exchange, not for domain decomposition.
