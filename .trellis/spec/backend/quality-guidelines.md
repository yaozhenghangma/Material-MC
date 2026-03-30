# Quality Guidelines

> Code quality standards for backend development.

---

## Overview

This project is a C++17 scientific simulation backend built with CMake and tested in CI by compilation.

Current quality gate in CI:

- `/Users/yma/Project/Material-MC/.github/workflows/main.yml`
  - installs MPI + build tools
  - runs `cmake .` and `make`

There is currently no configured formatter/linter/test framework in repo-level tooling.

---

## Forbidden Patterns

Based on existing architecture, avoid these patterns:

1. **Putting heavy logic in `src/MMC.cpp`**
   - Keep `MMC.cpp` as orchestration entrypoint.

2. **Mixing parsing/output code into compute kernels**
   - Keep parsing in `configure_in.cpp` / `structure_in.cpp`.
   - Keep output in `log.cpp` / `result_out.cpp` / `spin_out.cpp`.

3. **Breaking MPI root-only input contract**
   - In `src/MMC.cpp`, rank 0 parses input then broadcasts objects.

4. **Changing data structure fields without updating serialization**
   - `BaseSite`, `Lattice`, `MonteCarlo`, `Initialization` have cereal `serialize(...)` in `src/MC_structure.h`.

5. **Introducing frontend/web assumptions**
   - This repository has no frontend app code or HTTP service.

---

## Required Patterns

1. **Preserve header/source pairing and include guard style**
   - Examples: `log.h/.cpp`, `result_out.h/.cpp`.

2. **Register new source files in `CMakeLists.txt`**
   - Use `target_sources(MMC PRIVATE ...)`.

3. **Respect functional folder split**
   - Algorithms in `src/methods/`, extension logic in `custom/`.

4. **Keep runtime data in `Supercell`/`MonteCarlo` model**
   - Add shared state through `src/MC_structure.h` rather than ad hoc globals.

5. **Maintain MPI-safe flow in methods**
   - Follow existing gather/reorder patterns in classical/PTMC implementations.

---

## Testing Requirements (Current State)

Current required verification for changes:

1. **Build passes locally**
   - `cmake .`
   - `make`

2. **CI build remains green**
   - `.github/workflows/main.yml` build job

3. **For parser/output changes, run a sample simulation manually**
   - verify `log.txt`, `output.txt`, and optional `*.xsf` outputs are generated and readable.

Note: there is no formal unit-test suite yet.

---

## Code Review Checklist

- Does the change fit existing module boundaries (`src/`, `src/methods/`, `custom/`)?
- If structure fields changed, were cereal serialize lists updated in `src/MC_structure.h`?
- If new `.cpp`/`.h` files were added, were they added to `CMakeLists.txt`?
- Does MPI flow remain correct (root parse, broadcast, gather/reduce as needed)?
- Are outputs/logs still written via the existing output modules?
- Is behavior documented with comments consistent with current code style?

---

## Common Mistakes

- Forgetting to add new source files to CMake target sources.
- Editing hot-loop code without considering performance impact from extra logging/output.
- Adjusting TOML schema keys without matching parser updates in `ReadSettingFile`.
- Introducing partial feature logic in multiple files without clear single ownership.
