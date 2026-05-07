# Changelog

All notable changes to **Tadah.LAMMPS** (the `pair_tadah` LAMMPS pair-style integration) are documented here.
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

Nothing yet.

## [1.3.0-beta.1] — TBD

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

## [1.1.1] and earlier (condensed)

The `v1.0.0` → `v1.1.1` releases were the "internal playground" era. The LAMMPS pair-style integration matured during this period: CMake-driven LAMMPS build, integration tests, ghost-atom NN-list correctness fixes, and `D2_mJoin` / `compute_2b` wiring. `v1.1.1` was the last tag in this era and is the version on which downstream LAMMPS users have been running until this beta. No published changelog covers this era; treat anything before `1.3.0-beta.1` as undocumented in the modern sense.

> Note: Tadah.LAMMPS uses its own per-repo tag history (`v1.0.0`, `v1.1.1`); from `1.3.0-beta.1` onward it tracks the umbrella version.
