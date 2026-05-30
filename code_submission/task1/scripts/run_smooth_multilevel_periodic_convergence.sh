#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

export PROFILE=sinusoid
export PERIODIC_FLAGS="1 1"
export MAX_LEVEL=2
export STOP_TIME="${STOP_TIME:-0.05}"
export X0="${X0:-0.35}"
export Y0="${Y0:-0.40}"
export RHO0="${RHO0:-1.0}"
export AMP="${AMP:-0.05}"
export U0="${U0:-0.2}"
export V0="${V0:-0.15}"
export P0="${P0:-1.0}"
export SMOOTH_KX="${SMOOTH_KX:-1.0}"
export SMOOTH_KY="${SMOOTH_KY:-1.0}"
export PHIGRAD_SCALE="${PHIGRAD_SCALE:-0.128}"
export REGRID_INT="${REGRID_INT:-2}"
export N_ERROR_BUF="${N_ERROR_BUF:-2}"

exec "${SCRIPT_DIR}/run_smooth_amr_convergence.sh" "$@"
