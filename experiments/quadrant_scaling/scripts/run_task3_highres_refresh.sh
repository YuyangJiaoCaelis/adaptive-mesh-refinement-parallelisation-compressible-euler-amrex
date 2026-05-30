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

MPIEXEC="${MPIEXEC:-mpiexec}"
MPI_FLAGS="${MPI_FLAGS:---bind-to none --map-by :OVERSUBSCRIBE}"
read -r -a MPI_FLAG_ARR <<< "${MPI_FLAGS}"

UNIFORM_REPEATS="${UNIFORM_REPEATS:-3}"
AMR_REPEATS="${AMR_REPEATS:-3}"
UNIFORM_RANKS="${UNIFORM_RANKS:-1 2 4 8 16}"
AMR_RANKS="${AMR_RANKS:-1 2 4 8 16}"
RESULTS_DIR="${1:-${TASK_DIR}/results_task3_highres_refresh_$(date +%Y%m%d_%H%M%S)}"
mkdir -p "${RESULTS_DIR}"/{logs,raw,checks}

MASTER_LOG="${RESULTS_DIR}/checks/master.log"
CSV="${RESULTS_DIR}/checks/highres_refresh.csv"
SUMMARY_CSV="${RESULTS_DIR}/checks/highres_refresh_summary.csv"
SUMMARY_MD="${RESULTS_DIR}/checks/highres_refresh_summary.md"

echo "family,case_id,np,run_id,walltime_sec,level0_coverage_pct,level1_coverage_pct,level2_coverage_pct" > "${CSV}"

log() {
  echo "[INFO] $(date -Iseconds) $*" | tee -a "${MASTER_LOG}" >&2
}

parse_runtime() {
  local log_file="$1"
  grep -E "Run time =" "${log_file}" | tail -n 1 | awk '{print $4}'
}

parse_coverages() {
  local log_file="$1"
  awk '
    /Level[[:space:]]+[0-9]+[[:space:]]+[0-9]+ grids[[:space:]]+[0-9]+ cells[[:space:]]+[0-9.]+ % of domain/ {
      lev=$2
      pct=$(NF-3)
      cov[lev]=pct
    }
    END {
      l0=("0" in cov)?cov["0"]:"100.0"
      l1=("1" in cov)?cov["1"]:"0.0"
      l2=("2" in cov)?cov["2"]:"0.0"
      printf "%s,%s,%s\n", l0, l1, l2
    }
  ' "${log_file}"
}

should_keep_plotfile() {
  local family="$1"
  local np="$2"
  local run_id="$3"
  if [[ "${family}" == "uniform" && "${np}" -eq 1 && "${run_id}" == "r1" ]]; then
    return 0
  fi
  if [[ "${family}" == "amr" && "${np}" -eq 1 && "${run_id}" == "r1" ]]; then
    return 0
  fi
  return 1
}

run_uniform_case() {
  local np="$1"
  local run_id="$2"
  local case_id="uniform_n3072_np${np}"
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
    "amr.n_cell=3072 3072"
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
  echo "uniform,${case_id},${np},${run_id},${rt},,,">> "${CSV}"

  if ! should_keep_plotfile "uniform" "${np}" "${run_id}"; then
    find "${RESULTS_DIR}/raw" -maxdepth 1 -name "$(basename "${plot_prefix}")*" -exec rm -rf {} + >/dev/null 2>&1 || true
  fi
}

run_amr_case() {
  local np="$1"
  local run_id="$2"
  local case_id="amr_eff3072_np${np}"
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
    "adv.do_reflux=1"
    "geometry.is_periodic=0 0"
    "amr.max_level=2"
    "amr.n_cell=768 768"
    "amr.ref_ratio=2 2 2 2"
    "amr.regrid_int=2"
    "amr.n_error_buf=2"
    "amr.blocking_factor=16"
    "amr.max_grid_size=64"
    "amr.grid_eff=0.75"
    "amr.plot_int=1000000"
    "amr.plot_file=${plot_prefix}"
    "tagging.max_phierr_lev=-1"
    "tagging.max_phigrad_lev=2"
    "tagging.phigrad=0.08 0.04"
  )

  log "running ${case_id} ${run_id}"
  if [[ "${np}" -eq 1 ]]; then
    OMP_NUM_THREADS=1 "${solver_cmd[@]}" > "${log_file}" 2>&1
  else
    OMP_NUM_THREADS=1 "${MPIEXEC}" -n "${np}" "${MPI_FLAG_ARR[@]}" "${solver_cmd[@]}" > "${log_file}" 2>&1
  fi

  local rt cov l0 l1 l2
  rt="$(parse_runtime "${log_file}")"
  if [[ -z "${rt}" ]]; then
    echo "ERROR: runtime not found in ${log_file}" >&2
    exit 2
  fi
  cov="$(parse_coverages "${log_file}")"
  IFS=',' read -r l0 l1 l2 <<< "${cov}"
  echo "amr,${case_id},${np},${run_id},${rt},${l0},${l1},${l2}" >> "${CSV}"

  if ! should_keep_plotfile "amr" "${np}" "${run_id}"; then
    find "${RESULTS_DIR}/raw" -maxdepth 1 -name "$(basename "${plot_prefix}")*" -exec rm -rf {} + >/dev/null 2>&1 || true
  fi
}

for np in ${UNIFORM_RANKS}; do
  for rep in $(seq 1 "${UNIFORM_REPEATS}"); do
    run_uniform_case "${np}" "r${rep}"
  done
done

for np in ${AMR_RANKS}; do
  for rep in $(seq 1 "${AMR_REPEATS}"); do
    run_amr_case "${np}" "r${rep}"
  done
done

python3 - <<PY
import csv
import statistics
from collections import defaultdict
from pathlib import Path

csv_path = Path("${CSV}")
summary_csv = Path("${SUMMARY_CSV}")
summary_md = Path("${SUMMARY_MD}")

rows = list(csv.DictReader(csv_path.open()))
grouped = defaultdict(list)
coverages = defaultdict(list)
for row in rows:
    key = (row["family"], int(row["np"]))
    grouped[key].append(float(row["walltime_sec"]))
    if row["family"] == "amr":
        coverages[key].append(
            (
                float(row["level0_coverage_pct"] or 0.0),
                float(row["level1_coverage_pct"] or 0.0),
                float(row["level2_coverage_pct"] or 0.0),
            )
        )

uniform_p1_mean = statistics.mean(grouped[("uniform", 1)])

with summary_csv.open("w", newline="") as f:
    fieldnames = [
        "family", "np", "n", "mean_walltime_sec", "std_walltime_sec",
        "speedup_vs_uniform_p1_mean", "level0_coverage_pct",
        "level1_coverage_pct", "level2_coverage_pct"
    ]
    w = csv.DictWriter(f, fieldnames=fieldnames)
    w.writeheader()
    for key in sorted(grouped):
        family, np = key
        vals = grouped[key]
        row = {
            "family": family,
            "np": np,
            "n": len(vals),
            "mean_walltime_sec": statistics.mean(vals),
            "std_walltime_sec": statistics.stdev(vals) if len(vals) > 1 else 0.0,
            "speedup_vs_uniform_p1_mean": uniform_p1_mean / statistics.mean(vals),
            "level0_coverage_pct": "",
            "level1_coverage_pct": "",
            "level2_coverage_pct": "",
        }
        if family == "amr":
            l0 = statistics.mean(v[0] for v in coverages[key])
            l1 = statistics.mean(v[1] for v in coverages[key])
            l2 = statistics.mean(v[2] for v in coverages[key])
            row["level0_coverage_pct"] = l0
            row["level1_coverage_pct"] = l1
            row["level2_coverage_pct"] = l2
        w.writerow(row)

lines = []
lines.append("# Task 3 High-Resolution Refresh")
lines.append("")
lines.append("| Family | p | n | Mean runtime (s) | Std (s) | Speedup vs uniform p=1 mean |")
lines.append("|---|---:|---:|---:|---:|---:|")
for key in sorted(grouped):
    family, np = key
    vals = grouped[key]
    mean = statistics.mean(vals)
    std = statistics.stdev(vals) if len(vals) > 1 else 0.0
    speedup = uniform_p1_mean / mean
    lines.append(f"| {family} | {np} | {len(vals)} | {mean:.3f} | {std:.3f} | {speedup:.2f}x |")

summary_md.write_text("\\n".join(lines) + "\\n")
PY

log "results written to ${RESULTS_DIR}"
