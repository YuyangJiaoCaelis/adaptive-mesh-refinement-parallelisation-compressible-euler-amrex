#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
TASK_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
SOLVER_ROOT="${REPO_ROOT}/solver"
cd "${TASK_DIR}"

EXE="${SOLVER_ROOT}/build/main2d.gnu.MPI.ex"
if [[ ! -x "${EXE}" ]]; then
  echo "ERROR: executable not found: ${EXE}" >&2
  echo "Build first with: make -C \"${SOLVER_ROOT}/build\" -j8 AMREX_HOME=/path/to/amrex" >&2
  exit 1
fi

if [[ -z "${AMREX_HOME:-}" ]]; then
  echo "ERROR: set AMREX_HOME so fcompare.gnu.ex can be found." >&2
  exit 1
fi
FCOMPARE="${AMREX_HOME}/Tools/Plotfile/fcompare.gnu.ex"
if [[ ! -x "${FCOMPARE}" ]]; then
  echo "ERROR: fcompare not found: ${FCOMPARE}" >&2
  exit 1
fi

RESULTS_DIR="${1:-${TASK_DIR}/results_smooth_entropy_convergence_$(date +%Y%m%d_%H%M%S)}"
RAW_DIR="${RESULTS_DIR}/raw"
LOG_DIR="${RESULTS_DIR}/logs"
mkdir -p "${RAW_DIR}" "${LOG_DIR}"

STOP_TIME="${STOP_TIME:-0.2}"
X0="${X0:-0.35}"
Y0="${Y0:-0.40}"
RHO0="${RHO0:-1.0}"
AMP="${AMP:-0.2}"
U0="${U0:-0.2}"
V0="${V0:-0.15}"
P0="${P0:-1.0}"
SIGMA="${SIGMA:-0.08}"
X1=$(python3 - <<PY
x0=${X0}
u0=${U0}
t=${STOP_TIME}
print(f"{x0 + u0*t:.16g}")
PY
)
Y1=$(python3 - <<PY
y0=${Y0}
v0=${V0}
t=${STOP_TIME}
print(f"{y0 + v0*t:.16g}")
PY
)

latest_plotfile() {
  local prefix="$1"
  ls -d "${prefix}"* 2>/dev/null | grep -E "${prefix}[0-9]+$" | sort | tail -n 1 || true
}

for n in 128 256 512; do
  num_prefix="${RAW_DIR}/plt_smooth_num_n${n}"
  ref_prefix="${RAW_DIR}/plt_smooth_ref_n${n}"
  runlog="${LOG_DIR}/smooth_num_n${n}.log"
  reflog="${LOG_DIR}/smooth_ref_n${n}.log"
  rm -rf "${num_prefix}"* "${ref_prefix}"*

  OMP_NUM_THREADS=1 "${EXE}" inputs \
    "stop_time=${STOP_TIME}" \
    "prob.gamma=1.4" \
    "prob.ic=4" \
    "prob.print_ic=1" \
    "prob.x0=${X0}" \
    "prob.y0=${Y0}" \
    "prob.smooth_profile=1" \
    "prob.smooth_sigma=${SIGMA}" \
    "prob.smooth_rho0=${RHO0}" \
    "prob.smooth_amp=${AMP}" \
    "prob.smooth_u=${U0}" \
    "prob.smooth_v=${V0}" \
    "prob.smooth_p=${P0}" \
    "geometry.is_periodic=0 0" \
    "adv.first_order=0" \
    "adv.force_1d=0" \
    "adv.host_update=0" \
    "adv.do_reflux=0" \
    "adv.cfl=0.4" \
    "adv.v=0" \
    "amr.v=0" \
    "amr.max_level=0" \
    "amr.n_cell=${n} ${n}" \
    "amr.blocking_factor=8" \
    "amr.max_grid_size=${n}" \
    "amr.plot_int=1000000" \
    "amr.plot_file=${num_prefix}" \
    > "${runlog}" 2>&1

  OMP_NUM_THREADS=1 "${EXE}" inputs \
    "stop_time=0.0" \
    "prob.gamma=1.4" \
    "prob.ic=4" \
    "prob.print_ic=0" \
    "prob.x0=${X1}" \
    "prob.y0=${Y1}" \
    "prob.smooth_profile=1" \
    "prob.smooth_sigma=${SIGMA}" \
    "prob.smooth_rho0=${RHO0}" \
    "prob.smooth_amp=${AMP}" \
    "prob.smooth_u=${U0}" \
    "prob.smooth_v=${V0}" \
    "prob.smooth_p=${P0}" \
    "geometry.is_periodic=0 0" \
    "adv.first_order=0" \
    "adv.force_1d=0" \
    "adv.host_update=0" \
    "adv.do_reflux=0" \
    "adv.cfl=0.4" \
    "adv.v=0" \
    "amr.v=0" \
    "amr.max_level=0" \
    "amr.n_cell=${n} ${n}" \
    "amr.blocking_factor=8" \
    "amr.max_grid_size=${n}" \
    "amr.plot_int=1000000" \
    "amr.plot_file=${ref_prefix}" \
    > "${reflog}" 2>&1

  final_plot="$(latest_plotfile "${num_prefix}")"
  ref_plot="${ref_prefix}00000"
  if [[ -z "${final_plot}" || ! -d "${ref_plot}" ]]; then
    echo "ERROR: missing numerical or exact-reference plotfiles for N=${n}" >&2
    exit 2
  fi
  echo "DONE N=${n} final=${final_plot} ref=${ref_plot}"
done

python3 "${TASK_DIR}/analysis/smooth_entropy_convergence.py" \
  --results-dir "${RESULTS_DIR}" \
  --fcompare-exe "${FCOMPARE}"

echo "Smooth convergence results written to ${RESULTS_DIR}"
