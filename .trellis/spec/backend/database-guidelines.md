# Database Guidelines

> Data persistence conventions for this project.

---

## Overview

This repository currently has **no database layer**:

- no ORM
- no SQL/NoSQL client
- no migration system
- no persistence service module

The simulation is file-driven and in-memory:

- Input: `POSCAR` + `input.toml`
- Runtime state: `Supercell` / `MonteCarlo` objects in memory
- Output: `log.txt`, `output.txt`, `*.xsf`

Evidence in code:

- Input parsing: `/Users/yma/Project/Material-MC/src/configure_in.cpp`, `/Users/yma/Project/Material-MC/src/structure_in.cpp`
- Runtime data model: `/Users/yma/Project/Material-MC/src/MC_structure.h`
- Output writing: `/Users/yma/Project/Material-MC/src/result_out.cpp`, `/Users/yma/Project/Material-MC/src/spin_out.cpp`, `/Users/yma/Project/Material-MC/src/log.cpp`

---

## Query Patterns (Current Equivalent)

Since there is no DB, “query patterns” map to in-memory traversal patterns:

1. **Nested loop traversal over lattice/site containers**
   - Example: `Supercell::energy()` in `/Users/yma/Project/Material-MC/src/MC_structure.cpp`

2. **Neighbor adjacency pointer traversal during Hamiltonian evaluation**
   - Example: `Heisenberg*` functions in `/Users/yma/Project/Material-MC/src/Hamiltonian.cpp`

3. **Batch collection via MPI gather/reduction, then root aggregation**
   - Example: `/Users/yma/Project/Material-MC/src/methods/classical.cpp`
   - Example: `/Users/yma/Project/Material-MC/src/methods/parallel_tempering.cpp`

---

## Migrations

There is no schema migration tool in the current project.

Equivalent compatibility changes happen when editing:

- TOML key structure in `/Users/yma/Project/Material-MC/src/configure_in.cpp`
- POSCAR parsing assumptions in `/Users/yma/Project/Material-MC/src/structure_in.cpp`
- Runtime data structures in `/Users/yma/Project/Material-MC/src/MC_structure.h`

When changing any of these, update parsing and output logic consistently.

---

## Naming Conventions

Current “data-schema-like” naming follows existing C++/input conventions:

- TOML sections are PascalCase-style categories: `MonteCarlo`, `Lattice`, `Hamiltonian`, `Initialization`.
- TOML keys are snake_case: `start_temperature`, `temperature_points_number`, `exchange_step`.
- Runtime fields in C++ are snake_case: `temperature_step_number`, `neighbor_distance_square`, `super_exchange_parameter`.

Primary reference: `/Users/yma/Project/Material-MC/src/configure_in.cpp` and `/Users/yma/Project/Material-MC/src/MC_structure.h`.

---

## Anti-patterns / Common Mistakes

- Assuming transactional/DB semantics in this codebase (none exist).
- Changing input keys in `input.toml` without updating parser defaults in `ReadSettingFile`.
- Introducing hidden persistence side effects in compute kernels (`Hamiltonian.cpp`, `methods/*.cpp`).
- Treating output files as authoritative incremental state; outputs are run artifacts, not a database.
