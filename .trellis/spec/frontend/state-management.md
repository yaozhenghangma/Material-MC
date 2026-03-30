# State Management

> Frontend state management status in this project.

---

## Overview

Frontend state management is **not applicable** because no frontend application exists in this repository.

Current state is backend simulation state managed in C++ runtime objects.

---

## State Categories (Current Equivalent)

No frontend local/global/server state stores exist.

Equivalent backend runtime state lives in:

- `Supercell`, `Lattice`, `BaseSite`, `Site`, `MonteCarlo`, `Initialization`
- Defined in `/Users/yma/Project/Material-MC/src/MC_structure.h`

State initialization/parsing flow:

- CLI: `/Users/yma/Project/Material-MC/src/initialization.cpp`
- TOML: `/Users/yma/Project/Material-MC/src/configure_in.cpp`
- POSCAR: `/Users/yma/Project/Material-MC/src/structure_in.cpp`

---

## When to Use Global State

Not applicable to frontend.

Backend note: shared simulation data is passed through `Supercell`/`MonteCarlo` references and MPI communication rather than UI state containers.

---

## Server State

Not applicable currently.

No frontend server-state cache/query layer exists.

---

## Common Mistakes (Current Context)

- Describing Redux/Zustand/React Query patterns for a codebase that currently has none.
- Confusing backend simulation objects with frontend store abstractions.
