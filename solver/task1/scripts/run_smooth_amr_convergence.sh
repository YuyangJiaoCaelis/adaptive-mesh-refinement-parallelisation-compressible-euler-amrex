#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
TASK_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
CODE_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
cd "${TASK_DIR}"

EXE="${CODE_ROOT}/build/main2d.gnu.MPI.ex"
if [[ ! -x "${EXE}" ]]; then
  echo "ERROR: executable not found: ${EXE}" >&2
  exit 1
fi

RESULTS_DIR="${1:-${TASK_DIR}/results_smooth_amr_convergence_$(date +%Y%m%d_%H%M%S)}"
RAW_DIR="${RESULTS_DIR}/raw"
LOG_DIR="${RESULTS_DIR}/logs"
mkdir -p "${RAW_DIR}" "${LOG_DIR}"

STOP_TIME="${STOP_TIME:-0.05}"
X0="${X0:-0.35}"
Y0="${Y0:-0.40}"
RHO0="${RHO0:-1.0}"
AMP="${AMP:-0.05}"
U0="${U0:-0.2}"
V0="${V0:-0.15}"
P0="${P0:-1.0}"
SIGMA="${SIGMA:-0.12}"
PROFILE="${PROFILE:-gaussian}"
SMOOTH_KX="${SMOOTH_KX:-1.0}"
SMOOTH_KY="${SMOOTH_KY:-1.0}"
PHIGRAD_SCALE="${PHIGRAD_SCALE:-0.128}"
PHIGRAD_VALUES="${PHIGRAD_VALUES:-}"
REGRID_INT="${REGRID_INT:-2}"
N_ERROR_BUF="${N_ERROR_BUF:-2}"
MAX_LEVEL="${MAX_LEVEL:-1}"
PERIODIC_FLAGS="${PERIODIC_FLAGS:-0 0}"

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

for eff_n in 128 256 512; do
  base_n=$((eff_n / (2 ** MAX_LEVEL)))
  if [[ -n "${PHIGRAD_VALUES}" ]]; then
    phigrad=$(python3 - <<PY
values = [v.strip() for v in "${PHIGRAD_VALUES}".split(",") if v.strip()]
eff_n = ${eff_n}
mapping = {128: values[0], 256: values[1], 512: values[2]}
print(mapping[eff_n])
PY
)
  else
    phigrad=$(python3 - <<PY
scale=${PHIGRAD_SCALE}
base_n=${base_n}
print(f"{scale/base_n:.16g}")
PY
)
  fi
  max_grid_size=64
  if [[ "${base_n}" -lt "${max_grid_size}" ]]; then
    max_grid_size="${base_n}"
  fi

  num_prefix="${RAW_DIR}/plt_smooth_amr_eff${eff_n}"
  ref_prefix="${RAW_DIR}/plt_smooth_ref_eff${eff_n}"
  runlog="${LOG_DIR}/smooth_amr_eff${eff_n}.log"
  reflog="${LOG_DIR}/smooth_ref_eff${eff_n}.log"
  rm -rf "${num_prefix}"* "${ref_prefix}"*

  OMP_NUM_THREADS=1 "${EXE}" inputs \
    "stop_time=${STOP_TIME}" \
    "prob.gamma=1.4" \
    "prob.ic=4" \
    "prob.print_ic=0" \
    "prob.x0=${X0}" \
    "prob.y0=${Y0}" \
    "prob.smooth_profile=1" \
    "prob.smooth_sigma=${SIGMA}" \
    "prob.smooth_rho0=${RHO0}" \
    "prob.smooth_amp=${AMP}" \
    "prob.smooth_u=${U0}" \
    "prob.smooth_v=${V0}" \
    "prob.smooth_p=${P0}" \
    "prob.smooth_kx=${SMOOTH_KX}" \
    "prob.smooth_ky=${SMOOTH_KY}" \
    "prob.smooth_profile=$([[ \"${PROFILE}\" == \"gaussian\" ]] && echo 1 || echo 0)" \
    "geometry.is_periodic=${PERIODIC_FLAGS}" \
    "adv.first_order=0" \
    "adv.force_1d=0" \
    "adv.host_update=0" \
    "adv.do_reflux=1" \
    "adv.cfl=0.4" \
    "adv.v=0" \
    "amr.v=0" \
    "adv.diag_pressure_floor=1" \
    "amr.max_level=${MAX_LEVEL}" \
    "amr.n_cell=${base_n} ${base_n}" \
    "amr.ref_ratio=2 2 2 2" \
    "amr.regrid_int=${REGRID_INT}" \
    "amr.n_error_buf=${N_ERROR_BUF}" \
    "amr.blocking_factor=8" \
    "amr.max_grid_size=${max_grid_size}" \
    "amr.plot_int=1000000" \
    "amr.plot_file=${num_prefix}" \
    "tagging.max_phierr_lev=-1" \
    "tagging.max_phigrad_lev=${MAX_LEVEL}" \
    "tagging.phigrad=${phigrad}" \
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
    "prob.smooth_kx=${SMOOTH_KX}" \
    "prob.smooth_ky=${SMOOTH_KY}" \
    "prob.smooth_profile=$([[ \"${PROFILE}\" == \"gaussian\" ]] && echo 1 || echo 0)" \
    "geometry.is_periodic=${PERIODIC_FLAGS}" \
    "adv.first_order=0" \
    "adv.force_1d=0" \
    "adv.host_update=0" \
    "adv.do_reflux=0" \
    "adv.cfl=0.4" \
    "adv.v=0" \
    "amr.v=0" \
    "amr.max_level=0" \
    "amr.n_cell=${eff_n} ${eff_n}" \
    "amr.blocking_factor=8" \
    "amr.max_grid_size=${eff_n}" \
    "amr.plot_int=1000000" \
    "amr.plot_file=${ref_prefix}" \
    > "${reflog}" 2>&1

  echo "DONE effN=${eff_n}"
done

python3 "${SCRIPT_DIR}/smooth_amr_convergence.py" \
  --results-dir "${RESULTS_DIR}" \
  --stop-time "${STOP_TIME}" \
  --x0 "${X0}" \
  --y0 "${Y0}" \
  --rho0 "${RHO0}" \
  --amp "${AMP}" \
  --u0 "${U0}" \
  --v0 "${V0}" \
  --p0 "${P0}" \
  --sigma "${SIGMA}" \
  --profile "${PROFILE}" \
  --kx "${SMOOTH_KX}" \
  --ky "${SMOOTH_KY}" \
  --max-level "${MAX_LEVEL}"

echo "Smooth AMR convergence results written to ${RESULTS_DIR}"
