# Backend Development Guidelines

> Backend conventions for this repository (Material-MC).

---

## Overview

This project is a C++17 + MPI simulation backend (no web API/server layer).
The guides in this directory document **current, observed** conventions from the codebase so agents can match existing patterns.

---

## Guidelines Index

| Guide | Description | Status |
|-------|-------------|--------|
| [Directory Structure](./directory-structure.md) | Module organization and file layout | Done |
| [Database Guidelines](./database-guidelines.md) | Data persistence model (no DB/ORM in current repo) | Done |
| [Error Handling](./error-handling.md) | Error types, handling strategies | Done |
| [Quality Guidelines](./quality-guidelines.md) | Code standards, forbidden patterns | Done |
| [Logging Guidelines](./logging-guidelines.md) | Logging outputs and usage patterns | Done |

---

## Scope Note

These docs describe the current backend architecture centered around:

- `src/MMC.cpp` (program entry, MPI orchestration)
- `src/*.cpp` + `src/methods/*.cpp` (simulation logic)
- `custom/custom.cpp` (custom Hamiltonian extension point)
- `CMakeLists.txt` + `.github/workflows/main.yml` (build/CI expectations)

---

## Evidence Examples

- MPI orchestration and rank-0 parse/broadcast flow: `/Users/yma/Project/Material-MC/src/MMC.cpp`
- TOML input parsing and defaults: `/Users/yma/Project/Material-MC/src/configure_in.cpp`
- File-based logging/output conventions: `/Users/yma/Project/Material-MC/src/log.cpp`, `/Users/yma/Project/Material-MC/src/result_out.cpp`

---

## Common Mistakes

- Describing this repository as a web backend with routes/controllers.
- Adding generic backend rules that are not traceable to current files.
- Updating simulation architecture without syncing these guideline docs.

---

**Language**: All documentation is written in **English**.
