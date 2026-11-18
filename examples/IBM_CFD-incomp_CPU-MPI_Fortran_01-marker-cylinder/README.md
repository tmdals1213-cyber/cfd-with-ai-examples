# Incompressible marker IBM

Five-step MPI startup was completed; no developed wake or force accuracy claim.

```bash
python3 scripts/case.py build 01_ib_tester --jobs 8
python3 scripts/case.py run 01_ib_tester --input smoke
# Longer demonstration, explicitly selected:
python3 scripts/case.py run 01_ib_tester --input demo
```

MPI partition: `2 x 2 x 1` (4 ranks). `input.smoke` is deliberately short. A new machine's smoke output must be checked before longer runs. The publication wrapper/cold-start inputs are new portability work; historical results refer to the recorded original runs, not a fresh execution of these wrappers.

`src/` preserves source bytes from `ib_tester`. The input duration, output periods and checkpoint identity were adapted for a clean start. `provenance/source_manifest.json` binds copied source files. See [model limits](../../docs/MODELS_AND_LIMITS.md).
