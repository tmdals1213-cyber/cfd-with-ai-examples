# What changed from upstream examples

- Added EOS-based liquid-water enthalpy, gravity and grid-aligned solid conduction. A common interface conductance transfers equal/opposite heat; native low-Mach momentum/projection remains the base.
- Added separate liquid/vapor thermal fields, equilibrium interface heat flux and phase-change source coupling. The planar Stefan reference isolates this from general CHT/nucleation claims.
- Added prescribed thermally paid wall birth and diagnostics of gas centroid/contact area. This is an explicit model input, not a resolved microscopic nucleation law.
- Repaired a local open-boundary momentum/pressure consistency path in the 2D heat-driven study. Keep that code with its matching pressure/projection and boundary inputs; do not transplant only a velocity extrapolation line.
- Replaced a source-only phase inventory restriction with joint advection/source availability in the dielectric development. Effective mass transfer, expansion and latent/sensible accounting remain consistent. A disappearing phase transfers resolved sensible energy; accepted negative volumes and velocity clipping are forbidden. Finite-step coupling/temporal accuracy remains limited.
- Extended the thermal, momentum, phase transport, boundary and output paths into XYZ. Added true quadrilateral PLIC face clipping, three-direction owned staggered update ranges and 3D pressure gauge ownership. Native HYPRE offset indexing is `(stencil, coordinate)`.
- Added isolated portable builds/runs, source hashes, MPI-abort detection and short physical-history audits. This release changes packaging, durations/output periods and cold-start identities; it does not overwrite the historical solver states.

All custom case source is explicit in each directory. The recorded source hashes distinguish old example paths from the current 3D path. Do not apply the whole release as a patch to an unrelated solver.

2026-11-18 portfolio expansion:21entries covering vehicles, runners, cooling and explicit stopped/reference work. Added mixed-license root mapping, local GitHub validation, actual25–35ms preview and continuous-field renderer. Current3D v6 retains the separately tested output/moment repairs and raises only the thermal iteration budget; original residual criteria remain. Actual100ms completion is separate.
