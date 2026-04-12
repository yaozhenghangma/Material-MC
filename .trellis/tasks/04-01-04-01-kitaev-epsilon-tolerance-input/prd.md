# Task: 04-01-kitaev-epsilon-tolerance-input

## Overview
Make Kitaev geometry thresholds configurable through `input.toml` by introducing epsilon/tolerance input parameters with stable defaults, so KH bond-direction classification can be tuned without recompilation.

## Requirements
- Add `KitaevEpsilon` and `KitaevTolerance` (or project-consistent key names for KH epsilon/tolerance) to `input.toml` parsing for KH runs.
- Preserve backward compatibility by applying defaults when keys are omitted.
- Store parsed values in runtime data structures used by KH initialization logic.
- Replace hardcoded KH threshold constants in initialization/classification code with configured runtime values.
- Keep parser behavior consistent with existing `ReadSettingFile` patterns (`value_or(default)`, numeric validation/error style).
- Update KH example input and user-facing input schema documentation.

## Acceptance Criteria
- [ ] KH runs accept configurable epsilon/tolerance keys from `input.toml`.
- [ ] Missing keys fall back to documented default values and behavior remains unchanged from current defaults.
- [ ] KH classification in `src/initialization.cpp` uses configured runtime values instead of hardcoded constants.
- [ ] Any newly added runtime fields are correctly wired in shared structures/serialization as needed.
- [ ] `example/kh_minimal/input.toml` demonstrates the new keys (or clearly reflects default usage).
- [ ] `docs/source/input_toml_tutorial.rst` documents key names, defaults, and usage scope.

## Technical Notes
- Existing KH thresholds are currently hardcoded in `src/initialization.cpp`; implementation should route these through parsed configuration.
- Parser updates should be implemented in `src/configure_in.cpp` with established defaulting conventions.
- Runtime storage likely belongs in shared model structs in `src/MC_structure.h`; if new fields are added, keep cereal serialization aligned.
- Relevant KH debug paths (e.g., ambiguous/invalid direction grouping) are sensitive to epsilon/tolerance and should retain clear error diagnostics.

## Out of Scope
- Changes to non-KH model behavior or schemas.
- New KH coupling models (J/K/G/Gp semantics), per-bond overrides, or Hamiltonian feature expansion.
- Visualization-only or unrelated include/style refactors.
