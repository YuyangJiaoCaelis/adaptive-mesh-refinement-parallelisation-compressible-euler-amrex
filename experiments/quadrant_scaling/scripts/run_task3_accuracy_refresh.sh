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
  exit 1
fi

LOWERRES_DIR="${1:?usage: run_task3_accuracy_refresh.sh LOWERRES_REFRESH_DIR HIGHRES_REFRESH_DIR [RESULTS_DIR]}"
HIGHRES_DIR="${2:?usage: run_task3_accuracy_refresh.sh LOWERRES_REFRESH_DIR HIGHRES_REFRESH_DIR [RESULTS_DIR]}"
RESULTS_DIR="${3:-${TASK_DIR}/results_task3_accuracy_refresh_$(date +%Y%m%d_%H%M%S)}"

mkdir -p "${RESULTS_DIR}"/{logs,raw,checks}

AMR1536_PREFIX="${RESULTS_DIR}/raw/plt_amr_eff1536_np1_r1"
AMR1536_LOG="${RESULTS_DIR}/logs/amr_eff1536_np1_r1.log"

OMP_NUM_THREADS=1 "${EXE}" inputs \
  "prob.ic=3" \
  "prob.x0=0.5" \
  "prob.y0=0.5" \
  "stop_time=0.3" \
  "adv.first_order=0" \
  "adv.force_1d=0" \
  "adv.host_update=0" \
  "adv.debug_state=0" \
  "adv.v=0" \
  "amr.v=0" \
  "adv.do_reflux=1" \
  "geometry.is_periodic=0 0" \
  "amr.max_level=1" \
  "amr.n_cell=768 768" \
  "amr.ref_ratio=2 2 2 2" \
  "amr.regrid_int=2" \
  "amr.n_error_buf=2" \
  "amr.blocking_factor=16" \
  "amr.max_grid_size=64" \
  "amr.grid_eff=0.75" \
  "amr.plot_int=1000000" \
  "amr.plot_file=${AMR1536_PREFIX}" \
  "tagging.max_phierr_lev=-1" \
  "tagging.max_phigrad_lev=1" \
  "tagging.phigrad=0.08" \
  > "${AMR1536_LOG}" 2>&1

python3 "${TASK_DIR}/analysis/task3_fullfield_accuracy.py" \
  --uniform768-plot "$(find "${LOWERRES_DIR}/raw" -maxdepth 1 -type d -name 'plt_uniform_n768_np1_r1*' | sort | tail -n 1)" \
  --uniform1536-plot "$(find "${LOWERRES_DIR}/raw" -maxdepth 1 -type d -name 'plt_uniform_n1536_np1_r1*' | sort | tail -n 1)" \
  --amr1536-plot "$(find "${RESULTS_DIR}/raw" -maxdepth 1 -type d -name 'plt_amr_eff1536_np1_r1*' | sort | tail -n 1)" \
  --amr3072-plot "$(find "${HIGHRES_DIR}/raw" -maxdepth 1 -type d -name 'plt_amr_eff3072_np1_r1*' | sort | tail -n 1)" \
  --uniform3072-plot "$(find "${HIGHRES_DIR}/raw" -maxdepth 1 -type d -name 'plt_uniform_n3072_np1_r1*' | sort | tail -n 1)" \
  --out-dir "${RESULTS_DIR}/checks"

echo "${RESULTS_DIR}"
