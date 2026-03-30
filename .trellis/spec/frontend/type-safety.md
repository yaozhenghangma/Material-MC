# Type Safety

> Frontend type-safety status in this project.

---

## Overview

Frontend TypeScript type-safety guidance is **not applicable** in current repository state.

There is no TypeScript frontend code (`tsconfig.json`, `*.ts`, `*.tsx` absent at project level).

The project uses C++ type safety and explicit data structures.

---

## Type Organization (Current Equivalent)

Core simulation types are defined in C++ headers:

- `/Users/yma/Project/Material-MC/src/MC_structure.h` (classes and enums)
- `/Users/yma/Project/Material-MC/src/constants.h` (numeric constants)
- `/Users/yma/Project/Material-MC/src/Hamiltonian.h` (function interfaces)

Build enforces C++17:

- `/Users/yma/Project/Material-MC/CMakeLists.txt` (`set(CMAKE_CXX_STANDARD 17)`)

---

## Validation (Current Equivalent)

No runtime validation library like Zod/Yup is used.

Current validation style is parser-level defaults and limited parse exception handling:

- TOML defaults via `.value_or(...)` in `/Users/yma/Project/Material-MC/src/configure_in.cpp`
- parse exception catch in same file

---

## Common Patterns

- Explicit structs/classes for domain state (`Supercell`, `Site`, `MonteCarlo`).
- Enum-based mode selection (`HamiltonianType`, `ModelType`, `Methods`).
- Function signatures with references for shared mutable state.

---

## Forbidden Patterns (Frontend Context)

- Adding speculative TypeScript conventions before TypeScript frontend code exists.
- Treating this backend-only repository as if it already had TS app architecture.
