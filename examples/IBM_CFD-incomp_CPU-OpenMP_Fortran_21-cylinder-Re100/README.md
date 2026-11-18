# Re100 cylinder baseline

Re100 interpolated IBM, 160×96×1 cells. Historical20000steps to time200, maxCFL.45961, finite, no recorded Poisson nonconvergence. Not a new grid/time convergence claim.

The IBM forcing loops retained the FIRSTPRIVATE wall-velocity fix. Source and small historical setup are included; native field/checkpoint/preprocessing binaries are excluded. Build with `python3 scripts/legacy_case.py build 21_cylinder`; prepare with `python3 scripts/legacy_case.py prepare 21_cylinder --variant default`. Recreate matching IBM preprocessing before running.
