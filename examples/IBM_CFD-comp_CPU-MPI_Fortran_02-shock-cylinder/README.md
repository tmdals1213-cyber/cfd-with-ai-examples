# Compressible shock / ghost-point IBM

Five-step MPI startup was completed; no developed shock-interaction validation.

```bash
python3 scripts/case.py build 02_ib_shock --jobs 8
python3 scripts/case.py run 02_ib_shock --input smoke
# Longer demonstration, explicitly selected:
python3 scripts/case.py run 02_ib_shock --input demo
```

MPI partition: `2 x 2 x 1` (4 ranks). `input.smoke` is deliberately short. A new machine's smoke output must be checked before longer runs. The publication wrapper/cold-start inputs are new portability work; historical results refer to the recorded original runs, not a fresh execution of these wrappers.

`src/` preserves source bytes from `ib_shock`. The input duration, output periods and checkpoint identity were adapted for a clean start. `provenance/source_manifest.json` binds copied source files. See [model limits](../../docs/MODELS_AND_LIMITS.md).
