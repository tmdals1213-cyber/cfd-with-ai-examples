# Isolated runner vs inverted-eight pacemakers

The final campaign used a C++ OpenMP penalty/porosity IBM with explicit Vreman SGS. It is distinct from the earlier Fortran ghost-point multibody experiment.

```bash
python3 scripts/legacy_case.py build 15_runner_final --compiler g++
```

This builds the reusable kernel only; it does not start a runner calculation. A small non-Windows export-visibility preamble was added; original and release hashes are recorded. `src/` includes the kinematics and reconstruction/preconditioner mathematics. `reference_pipeline/` retains the final source orchestration for an agent to adapt; its exact-artifact manifests and external geometry/motion assets are intentionally not shipped. It must not be treated as an automatic fresh-run command.

`case.json` defines the actual formation, speed/cadence, domain and numerical conventions. `RECORDED_PRESSURE_SUMMARY.json` is a compact historical result. The final surface95 path excluded overlapping/internal union faces and fitted positive quadrature to constant/affine identities. Pressure sampling could still use solid-side cells, and viscous force was absent; total drag is not qualified. Do not publish the algebraic102.95% pressure-component change as physical drag reduction.

The old promotional videos and raw results were deleted at the user's earlier closeout request. The new short summary uses the recorded formation and pressure-component means; it is not a recovered time-resolved flow movie.

[Runner transfer prompt](../../prompts/RUNNERS.md) · [Pipeline setup/limits](../../docs/LEGACY_SETUP.md)
