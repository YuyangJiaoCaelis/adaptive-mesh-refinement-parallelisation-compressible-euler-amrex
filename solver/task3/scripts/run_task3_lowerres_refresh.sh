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

MPIEXEC="${MPIEXEC:-mpiexec}"
MPI_FLAGS="${MPI_FLAGS:---bind-to none --map-by :OVERSUBSCRIBE}"
read -r -a MPI_FLAG_ARR <<< "${MPI_FLAGS}"

RESULTS_DIR="${1:-${TASK_DIR}/results_task3_lowerres_refresh_$(date +%Y%m%d_%H%M%S)}"
mkdir -p "${RESULTS_DIR}"/{logs,raw,checks}
CSV="${RESULTS_DIR}/checks/lowerres_refresh.csv"
SUMMARY_MD="${RESULTS_DIR}/checks/lowerres_refresh_summary.md"
MASTER_LOG="${RESULTS_DIR}/checks/master.log"

echo "N,cores,run_id,walltime_sec" > "${CSV}"

log() {
  echo "[INFO] $(date -Iseconds) $*" | tee -a "${MASTER_LOG}" >&2
}

parse_runtime() {
  local log_file="$1"
  grep -E "Run time =" "${log_file}" | tail -n 1 | awk '{print $4}'
}

should_keep_plotfile() {
  local n="$1"
  local np="$2"
  local run_id="$3"
  [[ "${np}" -eq 1 && "${run_id}" == "r1" ]]
}

run_case() {
  local n="$1"
  local np="$2"
  local run_id="$3"
  local case_id="uniform_n${n}_np${np}"
  local log_file="${RESULTS_DIR}/logs/${case_id}_${run_id}.log"
  local plot_prefix="${RESULTS_DIR}/raw/plt_${case_id}_${run_id}"

  local -a solver_cmd=(
    "${EXE}" inputs
    "prob.ic=3"
    "prob.x0=0.5"
    "prob.y0=0.5"
    "stop_time=0.3"
    "adv.first_order=0"
    "adv.force_1d=0"
    "adv.host_update=0"
    "adv.debug_state=0"
    "adv.v=0"
    "amr.v=0"
    "adv.do_reflux=0"
    "geometry.is_periodic=0 0"
    "amr.max_level=0"
    "amr.n_cell=${n} ${n}"
    "amr.ref_ratio=2 2 2 2"
    "amr.regrid_int=2"
    "amr.n_error_buf=2"
    "amr.blocking_factor=16"
    "amr.max_grid_size=64"
    "amr.grid_eff=0.75"
    "amr.plot_int=1000000"
    "amr.plot_file=${plot_prefix}"
    "tagging.max_phierr_lev=-1"
    "tagging.max_phigrad_lev=0"
    "tagging.phigrad=0.08 0.04"
  )

  log "running ${case_id} ${run_id}"
  if [[ "${np}" -eq 1 ]]; then
    OMP_NUM_THREADS=1 "${solver_cmd[@]}" > "${log_file}" 2>&1
  else
    OMP_NUM_THREADS=1 "${MPIEXEC}" -n "${np}" "${MPI_FLAG_ARR[@]}" "${solver_cmd[@]}" > "${log_file}" 2>&1
  fi

  local rt
  rt="$(parse_runtime "${log_file}")"
  if [[ -z "${rt}" ]]; then
    echo "ERROR: runtime not found in ${log_file}" >&2
    exit 2
  fi
  echo "${n},${np},${run_id},${rt}" >> "${CSV}"

  if ! should_keep_plotfile "${n}" "${np}" "${run_id}"; then
    find "${RESULTS_DIR}/raw" -maxdepth 1 -name "$(basename "${plot_prefix}")*" -exec rm -rf {} + >/dev/null 2>&1 || true
  fi
}

for n in 768 1536; do
  for np in 1 2 4 8 16; do
    for rep in 1 2 3; do
      run_case "${n}" "${np}" "r${rep}"
    done
  done
done

python3 - <<PY
import csv
import statistics
from collections import defaultdict
from pathlib import Path

csv_path = Path("${CSV}")
rows = list(csv.DictReader(csv_path.open()))
grouped = defaultdict(list)
for row in rows:
    grouped[(int(row["N"]), int(row["cores"]))].append(float(row["walltime_sec"]))

lines = ["# Lower-Resolution Uniform MPI Refresh", ""]
lines.append("| N | cores | n | mean runtime (s) | std (s) |")
lines.append("|---:|---:|---:|---:|---:|")
for key in sorted(grouped):
    vals = grouped[key]
    lines.append(
        f"| {key[0]} | {key[1]} | {len(vals)} | {statistics.mean(vals):.6f} | "
        f"{statistics.stdev(vals) if len(vals) > 1 else 0.0:.6f} |"
    )
Path("${SUMMARY_MD}").write_text("\\n".join(lines) + "\\n")
PY

echo "${RESULTS_DIR}"
