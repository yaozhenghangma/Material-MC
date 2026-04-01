# Move misplaced includes from .cpp to matching .h

## Goal
Restore the project’s original include organization: declarations should carry their required header dependencies in the corresponding `.h` files, while `.cpp` files keep implementation-focused includes.

## Requirements
- Identify `.cpp` files that currently contain header includes that should belong to the matching module header.
- Move those includes into the corresponding `.h` files while preserving compile correctness.
- Keep include placement/style consistent with existing repository conventions.
- Do not change runtime behavior; this is an include-structure cleanup task only.
- After implementation is completed, update relevant `.trellis/spec/` guideline files to document/clarify the include-placement convention.

## Acceptance Criteria
- [ ] Misplaced includes are moved from affected `.cpp` files to their matching `.h` files.
- [ ] Affected `.h` files contain the dependencies needed by their declarations.
- [ ] Include order remains consistent with project style in touched files.
- [ ] Project build/check commands still pass after changes.
- [ ] Relevant spec file(s) are updated to reflect the include-placement rule.

## Technical Notes
- Scope is backend C++ code only.
- Prefer minimal edits: only move includes that violate the current project convention.
- Preserve public interface and behavior.
