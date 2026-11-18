# Project coverage

| Entry | Family | Scope |
|---|---|---|
| [01_ib_tester](../examples/IBM_CFD-incomp_CPU-MPI_Fortran_01-marker-cylinder) | NGA2 fundamentals | Five-step MPI startup was completed; no developed wake or force accuracy claim. |
| [02_ib_shock](../examples/IBM_CFD-comp_CPU-MPI_Fortran_02-shock-cylinder) | NGA2 fundamentals | Five-step MPI startup was completed; no developed shock-interaction validation. |
| [03_water_cht_2d](../examples/FVFD_CFD-lowMach-CHT_CPU-MPI_Fortran_03-water-natural-convection-2D) | thermal / phase change | Historical 10 s prototype; open outer reservoir, not the later closed-side tank. No phase change. |
| [04_stefan_2d](../examples/FVFD-VOF-PLIC_CFD-incomp-phasechange_CPU-MPI_Fortran_04-stefan-2D) | NGA2 fundamentals | One-grid analytic comparison completed to 0.5 s; not a general convergence proof. |
| [05_heated_bubble_2d](../examples/FVFD-VOF-PLIC_CFD-incomp-CHT-phasechange_CPU-MPI_Fortran_05-heated-bubble-2D) | thermal / phase change | Historical initially gas-free 100 W water case: wall birth, growth, detachment and rise observed through ~0.48 s; not a predictive boiling correlation. |
| [06_dielectric_pool_2d](../examples/FVFD-VOF-PLIC_CFD-incomp-CHT-phasechange_CPU-MPI_Fortran_06-dielectric-pool-2D) | thermal / phase change | Accepted evolution reached 0.1611 s, then stopped by user. Velocity temporal accuracy remains unresolved. |
| [07_immersion_3d](../examples/FVFD-VOF-PLIC_CFD-incomp-CHT-phasechange_CPU-MPI_Fortran_07-immersion-3D) | thermal / phase change | Current v6 retains the v5 moment repair and v4 output repair. Original25us failed-step crossing passed to37.2ms; native-CFL100us and thermal1000-iteration checks passed to39.6ms. Longer100ms continuation is separate; no full-cycle/temporal-accuracy qualification. |
| [08_ahmed](../examples/IBM_CFD-incomp_CPU-serial_Fortran_08-Ahmed-body) | vehicle | Re100 suspended coarse body; historical2000steps/time10; not a road-car accuracy benchmark. |
| [09_drivaer](../examples/IBM_CFD-incomp_CPU-OpenMP_Fortran_09-DrivAer-lowRe) | vehicle | Re500 coarse functional DrivAer test; final serial-to4thread continuation, no moving road. |
| [10_drivaer_les](../examples/IBM_CFD-incomp-LES_CPU-OpenMP_Fortran_10-DrivAer-wake) | vehicle | Re2000 Vreman wake demo; historical13141steps/time40; periodicity was not imposed. |
| [11_ev_pair](../examples/IBM_CFD-incomp-LES_CPU-OpenMP_Fortran_11-Model3-vs-IONIQ5) | vehicle | 100km/h public artist models, moving road/fixed wheels; qualitative comparison, not certified vehicle predictions. |
| [12_suv_pair](../examples/IBM_CFD-incomp-LES_CPU-OpenMP_Fortran_12-Urus-vs-GClass) | vehicle | 100km/h public artist models; current retained coarse direct-forcing route, moving-road approximation. |
| [13_car10](../examples/IBM_CFD-incomp-LES_CPU-OpenMP_Fortran_13-ten-vehicle-campaign) | vehicle | Ten artist-model bodies, common100km/h/time window; retain original force/operator semantics and shape limitations. |
| [14_turbine](../examples/IBM-moving_CFD-incomp-LES_CPU-OpenMP_Fortran_14-turbine) | other moving IBM | Prescribed rotor motion, not structural FSI. Historical startup/restart work; no new completed-rotor validation here. |
| [15_runner_final](../examples/IBM-penalty_CFD-incomp-LES_CPU-OpenMP_CPP-Python_15-solo-vs-eight-runners) | runner | Final solo and inverted8 matched2.501953125s demonstration; independent C++ penalty/Vreman path, distinct from LICA ghost-point interpolation. |
| [16_runner_four](../examples/IBM-ghost_CFD-incomp_CPU_Fortran-Python_16-four-runner-prototype) | runner | Historical four-person moving-IBM prototype. Reference modules only, not proof of a completed/stable final run. The final matched campaign used the separate C++ path. |
| [17_hull_stopped](../examples/FVM-VOF_CFD-incomp-FSI_CPU_CPP-Fortran_17-hull-stopped) | archived multiphysics | Stopped hull/VOF exploratory source. No full free-floating FSI completion or production force qualification. |
| [18_thermal_earlier](../examples/IBM_CHT_CPU_CPP-Fortran-Python_18-earlier-chip-thermal) | supporting historical work | Earlier LICA-era thermal prototypes, superseded by the current NGA2 thermal/phase-change cases. |
| [19_acoustic](../examples/FEM_CFD-acoustic_CPU_Fortran_19-SPECFEM-reference) | supporting historical work | Earlier source/transfer preparation; do not equate packaging with an acoustic solver execution. |
| [20_hairdryer](../examples/PREP_CFD-CHT_CPU_Python_20-hairdryer-preparation) | supporting historical work | Earlier handoff/inventory preparation; no completed hairdryer CFD/thermal outcome newly asserted. |
| [21_cylinder](../examples/IBM_CFD-incomp_CPU-OpenMP_Fortran_21-cylinder-Re100) | baseline verification | Re100 interpolated IBM, 160×96×1 cells. Historical20000steps to time200, maxCFL.45961, finite, no recorded Poisson nonconvergence. Not a new grid/time convergence claim. |

Earlier intermediate source variants are collapsed into representative working paths and explicit reference-only entries. The stopped hull/VOF work remains stopped. Acoustic/hairdryer entries record preparation rather than claiming a new solved case. Entire raw results and installed environments are excluded.
