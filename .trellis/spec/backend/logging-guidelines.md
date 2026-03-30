# Logging Guidelines

> How logging is done in this project.

---

## Overview

Logging is file-based and split by artifact purpose:

- `log.txt` for run metadata and progress
- `output.txt` for thermodynamic tabular results
- `*.xsf` for spin structure snapshots

Primary logging/output implementation files:

- `/Users/yma/Project/Material-MC/src/MMC.cpp`
- `/Users/yma/Project/Material-MC/src/log.cpp`
- `/Users/yma/Project/Material-MC/src/result_out.cpp`
- `/Users/yma/Project/Material-MC/src/spin_out.cpp`

The project uses `fmt::output_file` (from fmt) for file output.

---

## Log Levels (Current State)

There is no formal level enum (debug/info/warn/error).

Current practical separation:

1. **Progress/info messages**
   - `logger.print("Successfully ...")` in `src/MMC.cpp`

2. **Configuration summary detail**
   - `WriteLog(...)` and `WriteHamiltonian(...)` in `src/log.cpp`

3. **Error output**
   - startup parse errors go to `std::cerr` in `src/configure_in.cpp`

---

## Structured Logging

`log.txt` is human-readable plain text, not JSON.

Structure is section-based (examples in `src/log.cpp`):

- Lattice constants
- Monte Carlo method
- Model type
- Hamiltonian description
- Anisotropy and magnetic field
- Per-element cell information

Thermodynamic output is tabular with header in `output.txt` (`src/result_out.cpp`).

---

## What to Log

Current convention logs these events:

- Successful input processing and initialization (`src/MMC.cpp`)
- Selected algorithm/model/hamiltonian and lattice metadata (`src/log.cpp`)
- Final simulation completion and result output completion (`src/MMC.cpp`)

For new features, keep logs at similar granularity: setup summary + major phase transitions.

---

## What NOT to Log

- Per-spin/per-step verbose dumps inside Monte Carlo hot loops (performance and log size impact).
- Unbounded repeated logging in inner loops (`methods/*.cpp`).
- Sensitive data is not currently a concern in this scientific codebase, but avoid logging unrelated local environment details.

---

## Anti-patterns / Common Mistakes

- Printing detailed loop diagnostics each update step.
- Mixing user-facing progress text into result tables (`output.txt`).
- Writing key setup metadata only to stdout and not to `log.txt`.
- Creating inconsistent file formats across output writers.
