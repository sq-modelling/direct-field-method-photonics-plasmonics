#!/bin/sh
# SPDX-FileCopyrightText: 2026 Qiang Sun
# SPDX-License-Identifier: BSD-3-Clause
#
# Mesh-convergence study of manuscript Fig. 2(b): PEC sphere of radius a,
# k_out a = 4.5, icosahedral Q6 meshes of level L = 4, ..., 14, that is
# N = 40 L^2 + 2 = 642, ..., 7842 surface nodes.  The scattered field at
# (0, 0, 1.2a) is read from a 3 x 3 planar sample and compared with Mie theory
# by exercises/pec_sphere_mie.py.
#
# Usage, from the repository root after `make -C code/dfm_field_solver`:
#
#     sh exercises/run_sphere_convergence.sh [RUN_ROOT [LEVEL ...]]
#
# RUN_ROOT defaults to ../dfm-runs/sphere-convergence and the levels to
# 4 5 6 7 8 9 10 11 12 13 14.  Each level runs in a new directory
# RUN_ROOT/levelLL; an existing directory is never overwritten, because the
# solver appends to existing result files.  Set FIELD_SOLVER to use an
# executable built elsewhere (for example with BUILD_DIR).

set -eu

here=$(cd "$(dirname "$0")" && pwd)
repo=$(dirname "$here")
case_dir=$here/pec_sphere_convergence_ka4p5
exe=${FIELD_SOLVER:-$repo/code/dfm_field_solver/.build/field_solver}
run_root=${1:-$repo/../dfm-runs/sphere-convergence}
[ $# -gt 0 ] && shift
levels=${*:-4 5 6 7 8 9 10 11 12 13 14}

if [ ! -x "$exe" ]; then
    echo "field_solver not found at $exe; run 'make -C code/dfm_field_solver' first" >&2
    exit 1
fi
mkdir -p "$run_root"

for L in $levels; do
    d=$run_root/level$(printf '%02d' "$L")
    if [ -e "$d" ]; then
        echo "$d already exists; remove it or choose another RUN_ROOT" >&2
        exit 1
    fi
    mkdir "$d"
    cp "$case_dir/Input_Phys_EM.dat" "$case_dir/Input_Source_EM.dat" "$d/"
    # Set the mesh level in the 'Icshdrl' record of Input_Geom.dat.
    sed "s/^'Icshdrl', *[0-9][0-9]*,/'Icshdrl', $L,/" "$case_dir/Input_Geom.dat" > "$d/Input_Geom.dat"
    if ! grep -q "^'Icshdrl', $L," "$d/Input_Geom.dat"; then
        echo "could not set mesh level $L in $d/Input_Geom.dat" >&2
        exit 1
    fi
    echo "level $L: N = $((40 * L * L + 2)) nodes, running in $d"
    t0=$(date +%s)
    (cd "$d" && "$exe" > solver.log 2>&1)
    t1=$(date +%s)
    echo $((t1 - t0)) > "$d/elapsed_s.txt"
done

python3 "$here/pec_sphere_mie.py" point "$run_root"/level[0-9][0-9]
