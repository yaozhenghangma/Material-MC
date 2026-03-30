# Hook Guidelines

> Frontend hook conventions in this project.

---

## Overview

Not applicable currently: this repository has no frontend framework and no custom hook usage.

Evidence:

- C++ source tree under `/Users/yma/Project/Material-MC/src/`
- CMake build pipeline in `/Users/yma/Project/Material-MC/CMakeLists.txt`
- CI compiles C++ target only (`/Users/yma/Project/Material-MC/.github/workflows/main.yml`)

---

## Custom Hook Patterns

Not applicable currently.

No `use*` hook modules or frontend lifecycle/state hooks exist.

---

## Data Fetching

Not applicable currently.

There is no browser/client-side data fetching layer.
Simulation input is file-based (`input.toml`, `POSCAR`) parsed in backend code:

- `/Users/yma/Project/Material-MC/src/configure_in.cpp`
- `/Users/yma/Project/Material-MC/src/structure_in.cpp`

---

## Naming Conventions

Not applicable currently.

---

## Common Mistakes (Current Context)

- Assuming React hooks should be documented for this repository state.
- Mapping C++ utility functions to frontend hook patterns.
