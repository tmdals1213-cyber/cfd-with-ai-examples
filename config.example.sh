# Copy to config.local.sh and adapt paths. Execute from the repository root.
export NGA_HOME="$(pwd)/vendor/NGA2"
# Ubuntu multiarch packages need header/lib convenience prefixes:
# bash scripts/ubuntu_prefixes.sh
export HYPRE_DIR="$(pwd)/vendor/ubuntu-prefixes/hypre"
export LAPACK_DIR="$(pwd)/vendor/ubuntu-prefixes/lapack"
export FFTW_DIR="$(pwd)/vendor/ubuntu-prefixes/fftw"
export IRL_DIR="$HOME/cfd/deps/irl"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
# Open MPI on WSL1 only:
# export CFD_WSL1=1
