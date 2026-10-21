#!/usr/bin/env bash
# Clone and compile the CNES SOS radiative transfer code.
#
# The repository carries 3.5 GB of CKD absorption coefficients that this
# validation never reads (it runs without gas absorption), so the clone
# is shallow and lives outside the work tree, under SOS_ROOT.
#
# The vertical resolution is raised before compiling.  SOS ships with
# CTE_OS_NT_MIN at 100, and 101 layers are not enough when the
# molecule-to-aerosol ratio varies with altitude: its path reflectance
# is then 0.8 % off its own converged value.  SOS_LAYERS sets it; 500 is
# what leaves 0.02 %, and CTE_OS_NT, the array bound, allows 600.
#
# Usage:
#     pixi run build-sos
#     SOS_ROOT=/somewhere pixi run build-sos

set -euo pipefail

SOS_ROOT="${SOS_ROOT:-$HOME/.cache/adjeff/sos}"
SOS_LAYERS="${SOS_LAYERS:-500}"
REPO="https://github.com/CNES/RadiativeTransferCode-SOS"
EXE="$SOS_ROOT/exe/SOS_ABS_MAIN.exe"

# An existing build is reused only if it has the layer count asked for:
# a stock 100-layer SOS silently reads 0.8 % off on the path reflectance.
if [ -x "$EXE" ] && grep -q "^#define CTE_OS_NT_MIN $SOS_LAYERS\$" \
    "$SOS_ROOT/inc/SOS.h"; then
    echo "sos: $EXE is already built, with $SOS_LAYERS layers"
    exit 0
fi

if [ ! -d "$SOS_ROOT/src" ]; then
    echo "sos: cloning $REPO into $SOS_ROOT"
    mkdir -p "$(dirname "$SOS_ROOT")"
    git clone --depth 1 "$REPO" "$SOS_ROOT"
fi

command -v gfortran >/dev/null || {
    echo "sos: gfortran not found; run through 'pixi run build-sos'" >&2
    exit 1
}

echo "sos: setting CTE_OS_NT_MIN to $SOS_LAYERS layers"
sed -i -E "s/^#define CTE_OS_NT_MIN [0-9]+/#define CTE_OS_NT_MIN $SOS_LAYERS/" \
    "$SOS_ROOT/inc/SOS.h"
grep -q "^#define CTE_OS_NT_MIN $SOS_LAYERS\$" "$SOS_ROOT/inc/SOS.h" || {
    echo "sos: could not set CTE_OS_NT_MIN in inc/SOS.h" >&2
    exit 1
}

echo "sos: building with $(gfortran --version | head -1)"
# The Makefile reads SOS_ABS_ROOT and writes into obj/ and exe/ under it.
SOS_ABS_ROOT="$SOS_ROOT" make -C "$SOS_ROOT/gen" -f Makefile.gfortran

test -x "$EXE" || { echo "sos: build produced no executable" >&2; exit 1; }
echo "sos: built $EXE"
