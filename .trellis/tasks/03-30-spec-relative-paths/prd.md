# Normalize spec absolute paths to relative

## Goal
Replace absolute filesystem paths in Trellis spec documents (e.g. `/Users/yma/...`) with project-root relative paths so docs are portable across developers' machines.

## Requirements
- Identify spec documents containing absolute paths under `.trellis/spec/`.
- Replace absolute paths with repository-relative paths from project root.
- Preserve document meaning and references.
- Do not change unrelated guideline content.

## Acceptance Criteria
- [ ] No `/Users/yma` absolute path remains in `.trellis/spec/` docs.
- [ ] Path references in affected docs are relative to repository root.
- [ ] Only relevant spec docs are modified.

## Technical Notes
- Use existing path style in docs where possible.
- Keep markdown links and examples consistent after replacement.
