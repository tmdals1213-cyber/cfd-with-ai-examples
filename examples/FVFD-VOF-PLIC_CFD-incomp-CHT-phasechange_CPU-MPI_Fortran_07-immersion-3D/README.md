# True XYZ immersion-cooling prototype

Current v6 retains the v5 moment repair and v4 output repair. Original25us failed-step crossing passed to37.2ms; native-CFL100us and thermal1000-iteration checks passed to39.6ms. Longer100ms continuation is separate; no full-cycle/temporal-accuracy qualification.

```bash
python3 scripts/case.py build 07_immersion_3d --jobs 8
python3 scripts/case.py run 07_immersion_3d --input smoke
# A longer cold-start demonstration, explicitly selected:
python3 scripts/case.py run 07_immersion_3d --input demo
```

MPI partition:2x2x2, 80 cubed cells,20mm cube. The wrappers are portable cold starts; they do not import private checkpoints. Check a new machine's smoke output before longer runs.

`src/` contains v6: the v5 phase-moment repair, v4 output-array ownership repair and a thermal linear iteration budget of1000 instead of250. The thermal1e-11 and pressure1e-8 residual criteria are unchanged. The original point sink left a nonzero first moment while the gas volume nearly vanished, moving its normalized centroid outside the domain; its advective projection overflowed. The v5 depleted-phase source removes the phase-average specific first moment of the current advected phase. Positive source placement remains at the interface. Phase-volume rates and latent/sensible/pressure equations and checks are retained. This is an explicit first-order cell-integrated source closure, not an exact swept-interface Stefan construction. No centroid/velocity clipping or old-field fallback is used.

[Source patch](../../provenance/cpu3d_phase_moment_v5.patch) · [Build/kernel/CFD/output checks](../../provenance/PHASE_MOMENT_V5_CHECKS.json) · [Model limits](../../docs/MODELS_AND_LIMITS.md).

[20ms landscape introduction](../../media/IMMERSION_COOLING_LANDSCAPE_DRAFT01.mp4) · [Portrait](../../media/IMMERSION_COOLING_PORTRAIT_DRAFT01.mp4) · [Actual25–35ms flow preview](../../media/IMMERSION_FLOW_25_TO_35MS_PREVIEW.mp4). The introduction used historicalv3; the later preview usedv4. These movies are not relabelled asv5/100ms results.

[Thermal budget patch](../../provenance/cpu3d_thermal_iteration_v6.patch) · [V6 build/short-CFD checks](../../provenance/THERMAL_ITERATION_V6_CHECKS.json). The old thermal stop needed598iterations to meet the same criterion; the test reached39.6ms.
