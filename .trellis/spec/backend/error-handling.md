# Error Handling

> How errors are handled in this project.

---

## Overview

Error handling is currently **lightweight and imperative**:

- parsing code may use `try/catch` and terminate on failure
- CLI help path terminates early
- most compute/output functions return `int` status but usually always return success (`0` or `1`)
- there is no centralized error type hierarchy

Representative examples:

- TOML parse failure path: `/Users/yma/Project/Material-MC/src/configure_in.cpp`
- CLI help/exit path: `/Users/yma/Project/Material-MC/src/initialization.cpp`
- status-return-style functions: `/Users/yma/Project/Material-MC/src/spin_out.cpp`, `/Users/yma/Project/Material-MC/src/result_out.cpp`

---

## Error Types (Current State)

No custom error classes are defined in current code.

Observed handling types:

1. **Library exception catch**
   - `catch (const toml::parse_error& err)` in `ReadSettingFile`.

2. **Process termination for unrecoverable setup errors**
   - `exit(-1)` after TOML parse error.
   - `exit(0)` for usage/help output.

3. **Default/fallback behavior instead of explicit error**
   - unknown Monte Carlo method string defaults to classical (`ReadSettingFile`).

---

## Error Handling Patterns

### 1) Parse early on root rank, fail fast

`src/MMC.cpp` parses configuration on rank 0, then broadcasts data. Invalid TOML currently prints to stderr and exits.

### 2) Keep hot loops mostly exception-free

`methods/*.cpp`, `Hamiltonian.cpp`, and update kernels avoid exception handling in inner loops.

### 3) Use return codes for procedural steps

Many helper functions return `int` to indicate completion style, though callers rarely branch on return values.

---

## API Error Responses

Not applicable in this repository because there is no HTTP/API layer.

Equivalent “user-facing” error outputs are:

- `std::cerr` messages (e.g., TOML parse failure)
- console usage text for `-h`

---

## Anti-patterns / Common Mistakes

- Adding `exit()` in low-level reusable functions (prefer propagating status where feasible).
- Silently swallowing parse/IO errors without stderr/log message.
- Inconsistent default-fallback behavior for invalid input keys.
- Assuming exception-safe rollback semantics in MPI/compute code paths (none are implemented).

---

## Practical Guidance for New Changes

When adding new parsing or setup logic, follow current style:

1. Validate as early as possible (setup stage).
2. Emit clear stderr/log context (file/key/value where possible).
3. Avoid heavy exception use in inner Monte Carlo loops.
4. Keep behavior consistent with existing fail-fast startup model.
