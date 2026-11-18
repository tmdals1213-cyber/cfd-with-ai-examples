# Recorded results

These are historical observations from the original environment. New portability checks are reported separately in provenance/RELEASE_CHECKS.json.

- **IBM startup:** 4 MPI, 128×64×1, 5 steps to 0.1 case time; maximum CFL 0.447368. Startup evidence only.
- **Compressible IBM startup:** 4 MPI, 400×200×1, 5 steps to 0.01 nondimensional time; maximum CFL 0.358069. Positive density/pressure/temperature; no developed-shock validation.
- **Water CHT:** the historical 10 s transient prototype completed with finite/thermal/mass checks. Initial geometry and open reservoir differ from the later tank. No boiling.
- **Stefan:** 0.5 s, one grid, separate phase temperatures; layer-thickness relative discrepancy 6.108e-5. Detailed values are in provenance/recorded_results.json.
- **Heated vapor:** historical 100 W initially gas-free 2D study showed detachment between 0.400–0.420 s and about 8.56 mm tracked centroid rise to 0.480 s. Repackaged cold-start execution is a new run; do not borrow that result as proof of the portable input.
- **2D dielectric pool:** user-stopped at 0.1611179 s with a normal checkpoint near 0.1595179 s. Known velocity temporal accuracy limitation remains.
- **3D:** 512,000 cells, 8 MPI, 804 accepted steps to 20 ms, 41 min 31 s. Peak summed owned-rank RSS ≈6.85 GiB. Final energy-ledger error −3.65e-9 J, mass relative residual 1.61e-10, divergence-minus-source 1.44e-9 s⁻¹, pressure relative residual 8.09e-9. These are numerical execution checks, not predictive accuracy.

The 20→100 ms continuation remains separate from this frozen 20 ms showcase. Replace the video source manifest only after complete snapshots and their physical times are verified. Never silently relabel a 20 ms frame as a later result.

Intermediate output-bindingv4 repaired the saved-array references. The historical20ms source wasv3. V4 built successfully, restored the normal1000-step checkpoint with all8rank CSV files byte-equal, matched the pre-failure same-dt history exactly, crossed its output failure to25.2ms and exported temperature/latent_sink exactly asfloat32versions of actualdoubleCSV. Hardware traps and physical tolerances were retained. See [v4 audit](../provenance/OUTPUT_BINDING_V4_CHECKS.json). This does not establish a complete cooling cycle or guarantee the pending100ms run.

Latest execution update (2026-11-18): the original-dt v4 continuation again exited145 at trial1486 after last accepted t=0.037017788578745006s. All ranks were closed and normal step1100 was preserved. The earlier short output-reference repair is qualified only within its recorded scope; it does not resolve this remaining stop. An isolated same-condition native-fatal/stage replay is investigating it. No failed field is adopted, and no100ms completion is claimed.

Current source update: phase-momentv5 passed11 kernel checks, all8zero-step restore CSV byte comparisons, the original25us failure crossing to37.2ms, and exact temperature/latent/VOF/velocity export checks. The negative phase source now removes the current advected phase-average specific first moment; positive source placement and phase/heat/pressure rate equations remain. This is a first-order cell-integrated depletion closure, not exact swept-interface Stefan geometry or a temporal-accuracy/full-cycle result. See provenance/PHASE_MOMENT_V5_CHECKS.json and preserve historical failures. Check live ownership before any new run.

Latest source isv6. It retains the separately qualifiedv5 moment closure and extends only thermal linear maxiterations250 to1000. The original1e-11 criterion was met in598iterations; the same-operator short run reached39.6ms. The first accepted step matched priorv5 history exactly. See provenance/THERMAL_ITERATION_V6_CHECKS.json. The100ms continuation is a separate pending result.
