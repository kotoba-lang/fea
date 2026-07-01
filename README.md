# kotoba-lang/fea

Zero-dep portable `.cljc` — restored from the legacy `kami-engine/kami-cae`
Rust crate (deleted in kotoba-lang/kami-engine PR #82 "Remove Rust workspace
from kami-engine") as part of the **clj-wgsl migration** (ADR-2607010930,
`com-junkawasaki/root`).

KAMI CAE: FEA mesh generation, material library, boundary conditions,
solvers, and post-processing.

**Named `fea`, not `cae-solver`** — `kotoba-lang/cae-solver` already exists
as an unrelated repo (a `solve` multimethod contract for ROM/LBM CFD
backend dispatch, part of a different clean-sheet vehicle-design stack).
This repo restores the legacy `kami-cae` FEA engine under a collision-free
name.

| Namespace | Restored from | Purpose |
|---|---|---|
| `fea.mesh` | `mesh` | Nodes/elements (Beam2/Tri3/Tri6/Quad4/Tet4/Tet10/Hex8)/sets, mesh stats, regular box hex-mesher |
| `fea.material` | `material` | Material definitions (linear-elastic/hyperelastic/elasto-plastic) + built-in presets (steel/aluminum/titanium/concrete) |
| `fea.boundary` | `boundary` | DOF bitmasks + boundary condition variants (displacement/force/pressure/temperature/convection) |
| `fea.solver` | `solver` | Dense Cholesky linear-system solve + linear-static bar (Beam2) element assembly/BC-application/solve |
| `fea.postprocess` | `postprocess` | Result field ranges, point probing, color-map data export |

Depends on `kotoba-lang/engineer` for shared contracts.

## Status

Restored — all 5 modules ported from the original 939-line Rust `lib.rs`,
with all 4 original Rust unit tests mirrored 1:1 in `test/fea_test.cljc`
(+1 smoke test) — 6 tests / 20 assertions, 0 failures, including a full
1D-bar FEA solve (Cholesky decomposition through displacement/stress
computation) verified against the analytic solution to 1e-6 relative
error. Pure data + pure functions throughout; no IO/GPU. Matrices are
flat row-major vectors — coarse-grained per-analysis assembly/solve
(small dense matrices), not a per-frame hot loop.

`fea.solver`'s linear-static solver currently supports meshes composed
entirely of Beam2 (1-D bar) elements, matching the original's documented
scope (unit cross-section area; callers scale Young's modulus
accordingly).

## Develop

```bash
clojure -M:test
```
