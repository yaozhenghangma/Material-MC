# Task: migrate-scnlib-latest

## Overview
Upgrade vendored `scnlib` from the current legacy v0.4-based commit to the latest stable release, then adapt project usage so POSCAR parsing and build integration continue to work on Linux and macOS.

## Requirements
- Update the `submodules/scnlib` submodule pointer to the latest stable upstream release commit.
- Keep root build integration working (`CMakeLists.txt` still resolves/links expected scnlib target for `MMC`).
- Update project call sites/includes impacted by API or header layout changes (primarily `src/structure_in.cpp`, plus related headers).
- Preserve current parsing behavior for lattice vectors, labels/counts, and atom coordinates.
- Explicitly verify the parser call using positional format fields (`"{0} {1} {2}"`) remains valid after migration (or is replaced with equivalent behavior if upstream requires it).
- Validate build compatibility on Linux (CI path) and macOS (local build path).

## Acceptance Criteria
- [ ] `submodules/scnlib` points to latest stable release commit and repository metadata stays consistent (`.gitmodules` unchanged unless required).
- [ ] Project builds successfully with migrated scnlib using configured CMake integration.
- [ ] `src/structure_in.cpp` compiles cleanly with updated `scnlib` API and preserves existing parse semantics.
- [ ] Linux build path used by CI remains green (`.github/workflows/main.yml` expectations still satisfied).
- [ ] macOS build verification is executed and documented in task notes/output.
- [ ] No unrelated refactors are introduced outside migration-required changes.

## Technical Notes
- Current `scnlib` usage is concentrated in `src/structure_in.cpp` (limited API migration blast radius).
- `src/structure_in.h` and `src/configure_in.h` include `<scn/scn.h>` and must be kept compatible with upstream header structure.
- Root `CMakeLists.txt` currently links `scn::scn`; confirm target contract still matches upgraded submodule.
- Current CI verifies Linux; this task additionally requires explicit macOS verification.
- Local quality guidance prefers out-of-source CMake build while CI currently uses in-source style; validate at least one supported path end-to-end for each platform requirement.

## Out of Scope
- Rewriting parser architecture beyond required API adaptation.
- Introducing new input formats or changing POSCAR semantics.
- Non-scnlib dependency upgrades.
- Build-system redesign unrelated to keeping migration compatible.
