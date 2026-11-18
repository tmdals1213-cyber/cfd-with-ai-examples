# NGA2 subset notice

This notice describes only the seven NGA2 example subset. Other portfolio families and their separate terms are listed in NOTICE.md and docs/LICENSE_MAP.md.

# Attribution

NGA2: copyright (c) 2018-present Olivier Desjardins, MIT license (see licenses/MIT-NGA2.txt).
Base source: https://github.com/desjardi/NGA2 at `28a58d937af89dd506456048be1cad58f6d383a9`.

These examples contain NGA2 example/module code and project-specific modifications. The upstream copyright and permission notice are retained. Custom thermal, phase-inventory, 3D geometry, monitoring and presentation additions are distributed under the same MIT terms. This repository is an independent example collection; it is not an official NGA2 release.

IRL, HYPRE, Open MPI, LAPACK/BLAS, FFTW, CoolProp and Python dependencies retain their own licenses. Their source trees/binaries are not bundled. A small generated water-property table is bundled for the historical liquid-water example; generated with CoolProp 8.0.0 / HEOS::Water, 101325 Pa, 20–90 °C, 0.05 K spacing. Do not extrapolate it outside that model range.

No original LICA geometry, third-party CAD, SSH material, credentials, installed Linux image, binary or raw restart data is distributed here. Legacy runner notes and a generic pressure-integration utility are included separately from the NGA2 examples.
