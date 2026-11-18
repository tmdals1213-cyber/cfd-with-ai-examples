# 14-turbine

Prescribed rotor motion, not structural FSI. Historical startup/restart work; no new completed-rotor validation here.

Original LICA copyright/noncommercial notices are retained in every source header. The MIT license at the repository root does not relicense these files. See [license map](../../docs/LICENSE_MAP.md).

```bash
python3 scripts/legacy_case.py build 14_turbine --compiler gfortran --jobs 1
python3 scripts/legacy_case.py prepare 14_turbine --variant default
```

These commands build a new binary and copy only small setup files into ignored working directories. They do not launch a flow calculation. Supply the matching licensed geometry/SDF and IBM preprocessing artifacts separately; none are silently borrowed from an old field/restart. See [legacy setup](../../docs/LEGACY_SETUP.md).

`preprocess_reference/ibmpre_stl.f90` records the common 161×81×81 LES/no-interpolation preprocessor variant. Its compile-time dimensions/switches MUST be adapted to this case; it is not a universal ready-to-run preprocessor for every grid. `inputs/` preserves historical small setup files, which can request historical restart behavior; create and review an explicit cold-start configuration before a new run.

[Agent transfer](../../prompts/VEHICLES.md). Raw fields/checkpoints/binaries are excluded. No old remote or account controller is included.
