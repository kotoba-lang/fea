(ns fea
  "KAMI CAE — Computer-Aided Engineering: FEA mesh generation, boundary
  conditions, material library, solvers, and post-processing. Restored
  from the legacy kami-engine/kami-cae Rust crate (deleted in
  kotoba-lang/kami-engine PR #82 'Remove Rust workspace from kami-engine')
  as part of the clj-wgsl migration (ADR-2607010930, com-junkawasaki/root).

  Named `fea` (not `cae-solver`) to avoid collision with the pre-existing,
  unrelated kotoba-lang/cae-solver repo (a `solve` multimethod contract
  for ROM/LBM CFD backend dispatch — a different initiative entirely;
  discovered when this restoration was attempted under the 'cae-solver'
  name).

  One namespace per original Rust module:
    fea.mesh        — nodes/elements/sets, mesh stats, box mesher
    fea.material     — material definitions + built-in engineering presets
    fea.boundary     — DOF masks + boundary condition variants
    fea.solver       — dense Cholesky solve + linear-static bar assembly/solve
    fea.postprocess  — result field ranges, point probing, color-map export

  Zero-dep portable CLJC — pure data + pure functions, no IO/GPU. Coarse-
  grained per-analysis assembly/solve (small dense matrices), not a
  per-frame hot loop, matching ADR-2607010930's CLJC scope. Depends on
  kotoba-lang/engineer for shared contracts.")
