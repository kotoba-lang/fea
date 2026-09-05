# kotoba-fea

[![CI](https://github.com/kotoba-lang/fea/actions/workflows/ci.yml/badge.svg)](https://github.com/kotoba-lang/fea/actions/workflows/ci.yml)

**Finite-element analysis (FEA) domain logic in pure Clojure.** A
[kotoba-lang](https://github.com/kotoba-lang) capability library: mesh
generation, material library, boundary conditions, a dense-Cholesky
linear-static bar solver, and post-processing.

This is a from-scratch `.cljc` port of `kami-cae`, a Rust crate that lived
in `kami-engine` (`orgs/kotoba-lang/kami-engine`). `kami-cae` was deleted
from `kami-engine`'s working tree without being committed, as part of
retiring the Rust game-engine workspace in favor of pure-Clojure "kotoba"
authority repos (ADR-2607010000). This repo exists so that deletion loses
no domain knowledge — `kami-cae`'s full source is still recoverable via
`git show HEAD:kami-cae/src/lib.rs` in `kami-engine`'s history for as long
as that history exists, but the authority for this domain now lives here.

Named "fea" (not "cae"/"cae-solver") — `kotoba-lang/cae-solver` already
exists as an unrelated ROM/LBM CFD solve-dispatch contract.

No network, no I/O in any domain namespace (the one narrow exception,
`kotoba.fea.material`'s EDN-resource loader, is documented below). Pure
data + pure functions, portable `.cljc` across JVM / ClojureScript / SCI /
GraalVM.

## Maturity

| | |
|---|---|
| Role | capability |
| Tests | 44 tests, 98 assertions across 9 namespaces (JVM) and 35 tests, 67 assertions (portable nbb suite), all green |
| Solver scope | linear-static, `:beam2` (1-D bar, axial) + `:tet4` (3-D) elements; isotropic thermal strain from `:temperature` BCs + `:reference-temperature` |

## Thermal strain (tet4 + beam2)

`kotoba.fea.boundary/temperature` prescribes nodal temperatures; pass
`{:reference-temperature T0}` (Kelvin — caller-supplied, never assumed) as
the 4th argument to `solver/solve-linear-static`. The material model must
carry `:thermal-expansion` [1/K] (the built-in presets do). The solve
carries the isotropic initial strain eps_th = alpha (T − T0) as equivalent
nodal forces and subtracts it in stress recovery. tet4 results expose the
full 6-voigt mechanical stress via `:stress-voigt` — the von Mises scalar
alone reads 0 for the hydrostatic state that fully constrained thermal
expansion produces. `:pressure`/`:convection` BCs (constructed by
`kotoba.fea.boundary` but not implemented by this solver) are rejected
loudly (`:unsupported-bc-type`) instead of being silently dropped.

## What was ported
