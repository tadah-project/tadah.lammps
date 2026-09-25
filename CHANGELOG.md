# Changelog

All notable changes to **Tadah.LAMMPS** (the `pair_tadah` LAMMPS pair-style integration) are documented here.
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

### Fixed
- **`pair_coeff` now warns when a DM_REAM potential has no `EMBF0` line (BUG-08).** Such a potential is evaluated with the default `EMBF0 ZERO` (F(0) = 0 by construction), which moves every atomic energy of a potential written by tadah 1.5.0 or earlier by -F(0) (for a non-linear model, forces change too). Previously LAMMPS printed only the same neutral banner as for an explicit `EMBF0 ZERO`. The tadah notice meant for this case never showed up in LAMMPS: the tadah Logger defaults to ERROR and writes only to std::cout. Now every `pair_coeff` that loads such a potential prints a LAMMPS `WARNING`, on screen and in log.lammps, with the same text as the tadah notice, and the banner reads `EMBF0 ZERO (key absent: default ZERO, ...; add EMBF0 FREE to reproduce tadah <= 1.5.0 energies)`. Potentials that declare `EMBF0` print exactly what they printed before. Evaluation does not change. `LammpsTadah::embf0_declared` records whether the potential context as given stored the key, before DM_REAM materialises the default. **Rebuild required**: the `LammpsTadah` layout changed, and Tadah.HPO includes it.

## [1.5.1] — 2026-09-25

### Fixed
- **ML-TADAH now builds against LAMMPS with the runtime style registry (LAMMPS >= 2 Jun 2026).** LAMMPS commit `7ee0259e` ("style registry: stop generating the now-unused `style_*.h` and `packages_*.h`") moved pair/fix/compute styles from a compile-time `#include` list to a runtime registry: `styles/style_pair.h` is gone, and `GenerateStyleSources` now emits `styles/style_pair.cpp` which calls `Force::pair_styles().add_builtin(...)`. The integration hook appended `#include "pair_tadah.h"` to that header and hard-errored when it was missing, so **every configure against a LAMMPS newer than 2 Jun 2026 failed** with `ML-TADAH: Cannot find generated file: <build>/styles/style_pair.h`. Worse is the near miss: a LAMMPS configured *without* the hook links and starts perfectly, and only fails at run time with `ERROR: Unrecognized pair style 'tadah'` — which under `tadah hpo` surfaces as 100 % band-B "script failed" and no cause at all, because `LammpsSim` opens LAMMPS with `-screen none -log none`.
  The hook now calls LAMMPS's own `RegisterStyles()` on the `ML-TADAH` directory instead of patching a generated file. Both LAMMPS generations read the same global style properties, so one mechanism covers both, LAMMPS discovers the styles from the `PairStyle(...)` markers in our headers (no per-style bookkeeping here), and any style category we add later works without further changes. Verified end to end on **both** LAMMPS `2 Sep 2026` (generates `style_pair.cpp` with `add_builtin("tadah", &creator<PairTadah>)`) and `11 Feb 2026` (generates `style_pair.h` with `#include "pair_tadah.h"`): each configures, builds, lists `tadah` in `lmp -h`, and reaches `PairTadah::coeff()` on `pair_coeff * * EXTERN`.
  Timing is the subtle part and is commented in the hook: registration must happen after `include(StyleHeaderUtils)` defines `RegisterStyles` but before the generators consume the properties. `CMAKE_PROJECT_INCLUDE` runs hundreds of lines too early and `cmake_language(DEFER)` runs too late, so the hook watches `LAMMPS_STYLE_HEADERS_DIR`, which LAMMPS sets on the line immediately before those generators in both generations.
- **A configure-time assertion that the pair style was actually registered.** After LAMMPS finishes, the hook checks that a generated `style_pair.cpp`/`style_pair.h` really names `pair_tadah.h`, and hard-errors otherwise. A future LAMMPS reshuffle now stops the build instead of silently producing a binary that rejects `pair_style tadah` at run time.

### Changed
- `pair_coeff` prints the DM_REAM embedding reference (`EMBF0 ZERO` = F(0) = 0 by construction, the new default; `EMBF0 FREE` = raw 1.5.0 embedding). The EAMFAST kernel needs no change: F comes from `DM_Base::eam_embed`, which applies the reference. Rebuild required (DM_Base/DM_REAM layout). New DM_REAM collapse tests (Blip5 and negative-knot Knot5 fixtures): `eam_embed` equals w·calc_aed, F(0) == 0.0 under ZERO.
- `TADAH_LAMMPS_PAIR_STYLES` is **deprecated and ignored** — styles are discovered from the `PairStyle(...)` markers in the ML-TADAH headers. Setting it to a non-default value warns.

### Added
- `PairTadah::stage_extern(std::shared_ptr<tadah::LammpsTadah>)` / `clear_staged_extern()`: per-thread staging slot for the `EXTERN` potential handoff. Drivers that execute whole scripts through LAMMPS's `input->file()` (where per-line `load()` interleaving is impossible) stage the candidate potential once before the run; every `pair_coeff * * EXTERN ...` the script executes adopts it in `coeff()`. The slot lives outside the Pair object, so it survives `clear` cycles; precedence in `coeff()` is staged slot, then a pre-loaded `lt` (the 1.4.0 per-line driver pattern, unchanged), then a hard error.

### Fixed
- `pair_coeff * * EXTERN` with no staged or pre-loaded potential, and `pair_coeff` with an unloadable potential file, now raise a clean LAMMPS error. Previously `coeff()` discarded `load()`'s failure return and dereferenced the null potential — a segfault (this is how the legacy `input->file()`-based MPI worker died on EXTERN scripts).

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
