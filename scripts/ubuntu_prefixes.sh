#!/usr/bin/env bash
set -euo pipefail
here="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
arch="$(gcc -print-multiarch)"
base="$here/vendor/ubuntu-prefixes"
for name in hypre lapack fftw; do mkdir -p "$base/$name"; done
# These are local convenience symlinks, not a system install.
ln -sfn /usr/include/hypre "$base/hypre/include"
ln -sfn "/usr/lib/$arch" "$base/hypre/lib"
ln -sfn "/usr/lib/$arch" "$base/lapack/lib"
ln -sfn /usr/include "$base/fftw/include"
ln -sfn "/usr/lib/$arch" "$base/fftw/lib"
printf 'export HYPRE_DIR=%q\nexport LAPACK_DIR=%q\nexport FFTW_DIR=%q\n' "$base/hypre" "$base/lapack" "$base/fftw"
