#!/usr/bin/env bash
# Clone and compile the CNES SOS radiative transfer code.
#
# The repository carries 3.5 GB of CKD absorption coefficients that this
# validation never reads (it runs without gas absorption), so the clone
# is shallow and lives outside the work tree, under SOS_ROOT.
#
# Usage:
#     pixi run build-sos
#     SOS_ROOT=/somewhere pixi run build-sos

set -euo pipefail

SOS_ROOT="${SOS_ROOT:-$HOME/.cache/adjeff/sos}"
REPO="https://github.com/CNES/RadiativeTransferCode-SOS"
EXE="$SOS_ROOT/exe/SOS_ABS_MAIN.exe"

if [ -x "$EXE" ]; then
    echo "sos: $EXE is already built"
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

echo "sos: building with $(gfortran --version | head -1)"
# The Makefile reads SOS_ABS_ROOT and writes into obj/ and exe/ under it.
SOS_ABS_ROOT="$SOS_ROOT" make -C "$SOS_ROOT/gen" -f Makefile.gfortran

test -x "$EXE" || { echo "sos: build produced no executable" >&2; exit 1; }
echo "sos: built $EXE"
