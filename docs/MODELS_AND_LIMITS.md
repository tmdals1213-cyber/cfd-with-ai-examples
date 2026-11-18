# Models and retained limits

## Distinct numerical paths

The upstream marker-forcing incompressible IBM and ghost-point compressible IBM are separate examples. NGA2 is a Cartesian FV/FD research library; not every example uses IBM. VOF/PLIC tracks a fluid interface and is distinct from solid IBM. The chip cooling cases use grid-aligned solid conduction/flow masks rather than arbitrary production CPU CAD.

The water CHT case uses a fixed-thermodynamic-pressure EOS property table, conservative enthalpy and EOS-density buoyancy. It does not resolve acoustics, boiling or a water/air free surface. Its four open reservoir sides are the historical accepted prototype; the later tank-boundary development is not substituted silently.

Stefan is a fixed-temperature planar reference with separate phase temperatures. It is not a solid CHT or nucleation benchmark. The two-temperature result was 0.000705749609 m versus 0.000705706502 m analytic layer thickness at 0.5 s (relative thickness discrepancy 6.108e-5, growth-increment discrepancy 3.010e-4) on one grid only.

The heat-driven bubble case is a single-component liquid/vapor prototype with prescribed initial/wall birth, not spontaneous microscopic nucleation. Phase equilibrium uses fixed saturation properties. Local open-boundary momentum treatment was repaired; the original pressure/operator conventions are retained.

## Current dielectric / 3D model

FC-72-equivalent literature properties at 1 bar are fixed by phase; the selected source and complete property set are in [materials/fc72_proxy_1bar.json](../materials/fc72_proxy_1bar.json). This is not a validated commercial-fluid EOS or compressible two-phase thermodynamics. Liquid and vapor carry separate temperature/energy fields coupled to solid conduction, a saturation interface and geometric phase transport. Resolved ordinary transfer is bounded against actual post-advection phase availability and latent/sensible ledgers remain consistent. It can refuse difficult support/coupling rather than accept negative phase volume.

80³ case: 20 mm cube, 0.25 mm spacing, copper 8×8×1.5 mm heater, 3.79259 W, 12×12×1 mm cold plate at 45 °C, initial pool height 8 mm, initial top film 125 µm. Initial solid is 5 K above saturation (56.6 °C). Side/bottom walls are no-slip/adiabatic. Top outer vapor-reservoir ports balance phase conversion. A circular 1 mm prescribed wall-birth patch seeds a 0.5 µm equivalent embryo, paid by the energy ledger. This does not establish microscopic nucleation accuracy.

Top condensate may include reservoir-supplied vapor, and film motion includes an initial liquid film. A full sealed evaporation–condensation–return cycle has not been demonstrated. At 20 ms the heater has a thin attached vapor layer; detached bubble rise has not yet been established in 3D.

The source-only `phase_conversion_cfl` diagnostic may exceed 1 when advection replenishes a phase. It is not the hydrodynamic CFL. Do not remove the actual joint phase-volume guards or interpret the diagnostic alone as divergence.

## Known numerical issues and fixes

- Old source-only inventory rules stalled on vanishing phases. The current branch checks joint transport/source availability and consistently returns unused latent conversion energy to sensible energy. Velocity temporal accuracy is still unresolved in the 2D comparison.
- Native HYPRE stencil offsets are `(stencil, coordinate)`, and the pressure gauge must be owned in all three coordinates. Preserve the current 3D correction.
- Native staggered physical update/guard ranges exclude unused outer halo planes. Never turn valid active negative density into an accepted field; inspect actual owned/pressure-high faces.
- No detailed mechanism, long-time statistics, full grid/time convergence, heat-transfer correlation, leakage/drag or product-level CPU qualification is claimed.

## Legacy runner

The final independent C++ penalty/porosity IBM used Vreman SGS. It is separate from the early four-runner ghost-point experiment. Final surface95 quadrature excluded overlapping/interior union faces and fitted positive weights to closed-surface constant/affine identities. Its sampled pressure could still include solid-side cells, and no viscous traction was included. The final mean component was +8.3723N solo and −0.2466N paced over1.25–2.50s; it is not a total-drag reduction claim. Earlier raw full-mesh pressure sums are an older, weaker diagnostic. See the final runner case and prompts/RUNNERS.md.

Current source update: phase-momentv5 passed11 kernel checks, all8zero-step restore CSV byte comparisons, the original25us failure crossing to37.2ms, and exact temperature/latent/VOF/velocity export checks. The negative phase source now removes the current advected phase-average specific first moment; positive source placement and phase/heat/pressure rate equations remain. This is a first-order cell-integrated depletion closure, not exact swept-interface Stefan geometry or a temporal-accuracy/full-cycle result. See provenance/PHASE_MOMENT_V5_CHECKS.json and preserve historical failures. Check live ownership before any new run.

V6 changes only thermal linear iteration budget250 to1000; residual1e-11, pressure, physical operators and source closure are retained. Its598-iteration failure crossing reached39.6ms. Numerical execution checks are not whole-cycle or temporal-accuracy qualification.
