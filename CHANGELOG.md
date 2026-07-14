# Changelog

All notable changes to **Tadah.LAMMPS** (the `pair_tadah` LAMMPS pair-style integration) are documented here.
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

## [1.4.0] — 2026-07-14

### Added
- **EAM-collapse fast path** (`compute_2b_mb_half_eamfast`, new header `eam_collapse.h`): for a linear model whose descriptors are a radial two-body block plus an EAM-form many-body block, single species (the `tadah refit` family, e.g. `D2_Knot5` + `DM_REAM`), the potential collapses exactly to scalar pair/density/embedding functions. The pair style evaluates cubic-Hermite `phi`/`rho` tables sampled from the descriptors' own calc routines plus the analytic embedding, with scalar (1-double) ghost communication instead of descriptor-vector comm. ~7x faster force calls on the refit family; agrees with the generic kernel to ~1e-9 relative (a mid-node parity self-check is logged at load). Controlled by the `EAMFAST` config key, **default on**; set `EAMFAST false` in the potential/config for the bitwise generic path. Non-collapsible models (nonlinear, normalised, bias, multi-species, dimer, non-EAM-form MB) silently keep the generic kernel. Adjoint capture (`capture_adjoint_state`) with the fast path active is a hard error — such drivers must set `EAMFAST false`.
- `PairTadah::load(std::shared_ptr<tadah::LammpsTadah>)`: drivers (e.g. `tadah hpo`) can build ONE potential object per candidate and share it across `clear`/`pair_coeff EXTERN` cycles instead of reconstructing it per `pair_coeff`. The object carries per-simulation scratch, so sharing is valid only across sims that run sequentially (per-thread).

### Changed
- `PairTadah::lt` is now a `std::shared_ptr` (was `unique_ptr`); all existing `load` overloads keep their semantics.

## [1.3.0-beta.1] — 2026-05-20

A small but focused release: a per-species bias precomputation, the `i`/`j` swap fix, and the standard build / community / VERSION wiring.

### Added
- Version banner — `pair_tadah` now prints the Tadah! VERSION (read from the blueprint TOML at build time) from `PairTadah::coeff()` so that LAMMPS log files record exactly which version produced the run.
- `TADAH_ENABLE_ASAN` build option.

### Changed
- Per-species bias lookup precomputed at `initialize()` (instead of recomputed per call) — small but measurable per-step speed-up in LAMMPS runs that use the per-species intercept layout.

### Fixed
- Pair indexing: `i` / `j` swap corrected.

### Internal
- Community / governance pointer files added (umbrella canonical at `Tadah.MLIP/`).
- VERSION single-source-of-truth wired from blueprint TOML.
- CI: lint stage and tag-gated release-tarball stage added.
- Test build routed through `cmake-tools/cmake/TadahTesting.cmake` helpers (`tadah_declare_testing` / `tadah_add_tests`).

## [1.1.1] and earlier (condensed)

The `v1.0.0` → `v1.1.1` releases were the "internal playground" era. The LAMMPS pair-style integration matured during this period: CMake-driven LAMMPS build, integration tests, ghost-atom NN-list correctness fixes, and `D2_mJoin` / `compute_2b` wiring. `v1.1.1` was the last tag in this era and is the version on which downstream LAMMPS users have been running until this beta. No published changelog covers this era; treat anything before `1.3.0-beta.1` as undocumented in the modern sense.

> Note: Tadah.LAMMPS uses its own per-repo tag history (`v1.0.0`, `v1.1.1`); from `1.3.0-beta.1` onward it tracks the umbrella version.
