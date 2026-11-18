#!/usr/bin/env bash
set -euo pipefail
here="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
build="$here/build/07_immersion_3d/build/o/dp.gnu.opt.mpi"
casepath="$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["07_immersion_3d"]["directory"])' "$here/examples/catalog.json")"
work="$here/build/plic-face-checks"
mkdir -p "$work"
cd "$work"
gfortran -O0 -ffree-line-length-none -fcheck=all -ffpe-trap=invalid,zero,overflow -I"$build" "$here/$casepath/src/face_geometry3d.f90" "$here/tests/plic_faces.f90" "$build/precision_dp.o" -o check
./check
