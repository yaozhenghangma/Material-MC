# Frontend Development Guidelines

> Frontend status for this repository.

---

## Overview

This repository currently has **no application frontend** (no React/Vue/Svelte web UI, no TypeScript app, no component tree).

The codebase is backend/simulation oriented (C++ + MPI).

---

## Guidelines Index

| Guide | Description | Status |
|-------|-------------|--------|
| [Directory Structure](./directory-structure.md) | Frontend module layout | Not applicable (no frontend code) |
| [Component Guidelines](./component-guidelines.md) | Component patterns, props, composition | Not applicable (no frontend code) |
| [Hook Guidelines](./hook-guidelines.md) | Custom hooks, data fetching patterns | Not applicable (no frontend code) |
| [State Management](./state-management.md) | Local/global/server state | Not applicable (no frontend code) |
| [Quality Guidelines](./quality-guidelines.md) | Frontend code quality standards | Not applicable (no frontend code) |
| [Type Safety](./type-safety.md) | Type patterns, validation | Not applicable (no frontend code) |

---

## Evidence

- C++ backend entrypoint: `src/MMC.cpp`
- CMake-based build target: `CMakeLists.txt`
- CI builds C++ target only: `.github/workflows/main.yml`
- No frontend app/tooling files found (`package.json`, `tsconfig.json`, `*.tsx`, `*.jsx`)

---

## Future Note

If a frontend is added later, replace the “Not applicable” docs with real conventions from that codebase.

---

## Common Mistakes

- Copying generic frontend rules into this repository without any frontend files.
- Treating C++ simulation modules as if they were frontend components/hooks.
- Forgetting to update these docs after introducing actual frontend tooling.
