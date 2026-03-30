# Directory Structure

> How backend code is organized in this project.

---

## Overview

This repository is a **single C++ simulation backend** (MPI Monte Carlo), not a web backend with routes/controllers.
Code is organized by simulation responsibility: input parsing, model definitions, physics kernels, MC methods, and output.

---

## Directory Layout

```text
/Users/yma/Project/Material-MC/
├── CMakeLists.txt                  # Build target + source registration
├── src/
│   ├── MMC.cpp                     # Program entry and MPI orchestration
│   ├── MC_structure.h/.cpp         # Core runtime data model
│   ├── configure_in.h/.cpp         # input.toml parsing
│   ├── structure_in.h/.cpp         # POSCAR parsing
│   ├── initialization.h/.cpp       # CLI parsing + cell/neighbor initialization
│   ├── Hamiltonian.h/.cpp          # Built-in Hamiltonian kernels
│   ├── random_function.h/.cpp      # RNG helpers
│   ├── rotation.h/.cpp             # Spin rotation math
│   ├── log.h/.cpp                  # log.txt output
│   ├── result_out.h/.cpp           # output.txt thermodynamic table
│   ├── spin_out.h/.cpp             # XSF spin structure output
│   └── methods/
│       ├── local_update.h/.cpp     # Local Metropolis updates
│       ├── classical.h/.cpp        # Classical MC workflow
│       └── parallel_tempering.h/.cpp # PTMC workflow
├── custom/
│   └── custom.h/.cpp               # User custom Hamiltonian extension point
├── submodules/                     # Third-party libs (fmt, toml++, cereal, etc.)
└── .github/workflows/main.yml      # CI build (cmake + make)
```

---

## Module Organization

Current organization is functional and pipeline-driven:

1. **Entry + orchestration**
   - `src/MMC.cpp` handles MPI init, root-only input parsing, broadcast, method dispatch, and output.

2. **Input/parsing layer**
   - `src/configure_in.cpp` parses `input.toml`.
   - `src/structure_in.cpp` parses `POSCAR`.
   - `src/initialization.cpp` parses CLI options and builds runtime topology.

3. **Core model layer**
   - `src/MC_structure.h` defines `Supercell`, `Lattice`, `BaseSite`, `Site`, `MonteCarlo`, etc.

4. **Physics/update kernels**
   - `src/Hamiltonian.cpp` and `custom/custom.cpp` provide energy kernels.
   - `src/methods/local_update.cpp` provides local update rules.

5. **Sampling workflows**
   - `src/methods/classical.cpp` and `src/methods/parallel_tempering.cpp` implement algorithm-level loops and MPI data gathering.

6. **Output layer**
   - `src/log.cpp`, `src/result_out.cpp`, `src/spin_out.cpp` write run artifacts.

---

## Naming Conventions

Observed conventions in this repository:

- Files use **snake_case**: `configure_in.cpp`, `result_out.cpp`.
- Header/source pairs share the same base name: `log.h` + `log.cpp`.
- Monte Carlo algorithm implementations live under `src/methods/`.
- Public functions are commonly PascalCase/CamelCase: `ReadSettingFile`, `InitializeSupercell`, `ClassicalMonteCarlo`.
- Header guards are uppercase macros: `#ifndef CONFIGURE_IN`, `#ifndef PARALLEL_TEMPERING`.

---

## Anti-patterns / Common Mistakes

- Putting new algorithm logic directly in `src/MMC.cpp` instead of `src/methods/`.
- Adding parser logic into unrelated files (e.g., parsing in `Hamiltonian.cpp` instead of `configure_in.cpp` / `structure_in.cpp`).
- Writing new output formats inside compute kernels instead of `result_out.cpp` / `spin_out.cpp` / `log.cpp`.
- Bypassing `CMakeLists.txt` source registration when adding new `.cpp` files.

---

## Examples

- Program orchestration: `/Users/yma/Project/Material-MC/src/MMC.cpp`
- Input parsing split: `/Users/yma/Project/Material-MC/src/configure_in.cpp`, `/Users/yma/Project/Material-MC/src/structure_in.cpp`
- Method-specific algorithms: `/Users/yma/Project/Material-MC/src/methods/classical.cpp`, `/Users/yma/Project/Material-MC/src/methods/parallel_tempering.cpp`
- Extension point: `/Users/yma/Project/Material-MC/custom/custom.cpp`
