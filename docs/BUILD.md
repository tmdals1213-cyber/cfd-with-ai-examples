# Build once, select one case

Tested historical toolchain: Ubuntu 24.04 / WSL1, GNU Fortran 13.3, Open MPI 4.1.6, HYPRE 2.28, LAPACK 3.12. Use a native Linux filesystem for builds/runs. Python 3.10+ is required by the portable controller. This collection does not install a Linux distribution or alter BIOS settings.

```bash
sudo apt-get update
sudo apt-get install git make cmake gfortran g++ python3 libopenmpi-dev openmpi-bin libhypre-dev liblapack-dev libfftw3-dev libeigen3-dev
bash scripts/fetch_nga2.sh
bash scripts/ubuntu_prefixes.sh
cp config.example.sh config.local.sh
# Edit config.local.sh: NGA_HOME must be the fetched native Linux path.
source config.local.sh
python3 scripts/case.py list
python3 scripts/case.py build 01_ib_tester --jobs 8
python3 scripts/case.py run 01_ib_tester --input smoke
```

The Ubuntu prefix helper points HYPRE include to `/usr/include/hypre` and each library directory to the compiler multiarch library directory. It creates only local symlinks under ignored vendor/. If using a custom install, point to its prefix instead. Inspect the emitted link command on another distribution.

## IRL for VOF / phase-change examples

Historical IRL commit: `062aa7659a6a46e3ae13de6be8b3fa787902b8c6`.

```bash
git clone https://github.com/robert-chiodi/interface-reconstruction-library.git irl-source
git -C irl-source checkout --detach 062aa7659a6a46e3ae13de6be8b3fa787902b8c6
cmake -S irl-source -B irl-build -DCMAKE_BUILD_TYPE=Release -DIRL_BUILD_FORTRAN=ON -DBUILD_TESTING=OFF -DCMAKE_INSTALL_PREFIX="$HOME/cfd/deps/irl" -DCMAKE_CXX_COMPILER=g++ -DCMAKE_Fortran_COMPILER=gfortran
cmake --build irl-build --parallel 8
cmake --install irl-build
export IRL_DIR="$HOME/cfd/deps/irl"
```

NGA2 expects `include/irl_fortran` modules and `libirl_fortran`, `libirl_c`, `libirl`, `libabsl_all`. This IRL recipe is transcribed from the actual local setup, but a clean OS installation has not been repeated for this release. CMake may need network access to upstream dependencies.

## Execution

The wrapper copies the selected source into an isolated build directory and each run into a fresh output directory. It never overwrites the original examples. It refuses changed upstream commits, nonmatching source hashes and a second solver managed by the same checkout. It records source/input/binary identity and starts MPI at the partition size declared in the input; it does not choose a new decomposition silently.

`CFD_WSL1=1` enables the historical Open MPI shared-memory compatibility flags. Do not carry those MCA options to a different MPI implementation blindly. Set MPI binding with `CFD_MPI_EXTRA` if needed; default `--bind-to none` avoids assuming a particular CPU topology. Keep OpenBLAS/OMP at one thread per MPI rank.

No memory ceiling is imposed. The 80³ case used about 6.9 GiB summed rank RSS in the historical 20 ms calculation. That is a measured example, not a portable memory guarantee.

The runner checks logs for MPI_ABORT even if the OS return code is zero, and checks the custom history for finite values and recorded conservation/residual limits. Upstream startup cases need their `monitor/` values inspected. The runner's success is an execution check, not physics validation. External schedulers and other checkouts must be coordinated separately.

Short CLI IDs (such as `07_immersion_3d`) map to the descriptive folder names in examples/catalog.json. Each case keeps its method, physics, CPU-MPI resource and Fortran language visible in the folder name.

## Small checks without another CFD run

```bash
python3 -B scripts/verify_release.py
python3 -B tests/check_controller.py
# After building 07_immersion_3d:
bash scripts/check_plic_faces.sh
```

The kernel tests cover full/empty/half/oblique clipped faces and all three face orientations. Audit tests cover normal history, unfinished time, excess mass error and MPI_ABORT with an apparently successful OS exit.

The source manifest preserves original file-byte hashes. Git attributes disable automatic line-ending conversion; source files may retain historical LF or CRLF. The tested GNU toolchain accepts these source files. Shell setup scripts are authored with LF. Do not normalize copied source bytes without updating their release hashes deliberately.

Release preparation on 2026-11-18: all seven isolated case builds passed against the pinned NGA2 source and existing local dependencies. Six PLIC face-kernel checks and four controller audit tests passed. No second CFD campaign was launched for packaging. A new clean Linux installation and the new cold-start runs are still separate checks; see provenance/RELEASE_CHECKS.json.

The v5 depleted-phase moment kernel needs only GNU Fortran, not MPI/HYPRE/IRL: `python3 -B scripts/check_phase_moments.py`. Its11checks include actual failure operands and subnormal support. They are also included in the upload validation workflow. The prior seven-case build matrix records historicalv3 for case07; currentv5 build/short-CFD proof is provenance/PHASE_MOMENT_V5_CHECKS.json.
