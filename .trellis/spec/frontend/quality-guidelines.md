# Quality Guidelines

> Frontend quality standards status in this project.

---

## Overview

Frontend quality rules are currently **not applicable** because there is no frontend codebase.

Current CI quality gate is backend build verification:

- `.github/workflows/main.yml`
  - installs toolchain
  - runs CMake + make build steps

---

## Forbidden Patterns (Current Context)

- Adding frontend lint/test requirements that have no source files to enforce them.
- Documenting accessibility/UI review checks for a non-existent UI layer.
- Claiming React/TypeScript tooling exists when repository evidence shows C++/CMake only.

---

## Required Patterns (Current Context)

- Keep frontend guideline docs explicit about non-applicability until frontend code is introduced.
- Base documentation on observed repository evidence (files, build config, CI pipeline).

---

## Testing Requirements

Frontend testing requirements are not applicable currently.

Repository-level quality checks today are backend checks:

1. Local backend build succeeds from `build/` (out-of-source): `cmake .. && make`
2. GitHub Actions build succeeds (`.github/workflows/main.yml`)

---

## Code Review Checklist (Frontend Section)

For now, reviewers should verify:

- Frontend docs correctly state “no frontend code present”.
- Any future frontend guidance is only added after frontend files/tooling are introduced.

---

## Common Mistakes

- Copying generic frontend standards into this repo without evidence.
- Forgetting to update these “not applicable” sections when a real frontend is eventually added.
