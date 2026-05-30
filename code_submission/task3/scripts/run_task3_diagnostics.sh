#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
TASK_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
CODE_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
cd "${TASK_DIR}"

EXE="${CODE_ROOT}/build/main2d.gnu.MPI.ex"
MPIEXEC="${MPIEXEC:-mpiexec}"
MPI_FLAGS="${MPI_FLAGS:---bind-to none --map-by :OVERSUBSCRIBE}"
read -r -a MPI_FLAG_ARR <<< "${MPI_FLAGS}"

RESULTS_DIR="${1:-${TASK_DIR}/results_task3_diagnostics_$(date +%Y%m%d_%H%M%S)}"
mkdir -p "${RESULTS_DIR}"/{logs,raw,checks}

CSV="${RESULTS_DIR}/checks/task3_diagnostics.csv"
SUMMARY_MD="${RESULTS_DIR}/checks/task3_diagnostics.md"
echo "case_id,np,walltime_sec,fillpatch_sec,advance_sec,reflux_sec,tagging_sec,fillpatch_pct,advance_pct,reflux_pct,tagging_pct,pressure_floor_stage1,pressure_floor_stage2,pressure_floor_total" > "${CSV}"

parse_runtime() {
  local log_file="$1"
  grep -E "Run time =" "${log_file}" | tail -n 1 | awk '{print $4}'
}

run_case() {
  local case_id="$1"
  local np="$2"
  local log_file="${RESULTS_DIR}/logs/${case_id}.log"
  local plot_prefix="${RESULTS_DIR}/raw/plt_${case_id}"
  rm -rf "${plot_prefix}"*

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
    "adv.diag_timers=1"
    "adv.diag_pressure_floor=1"
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

  if [[ "${np}" -eq 1 ]]; then
    OMP_NUM_THREADS=1 "${solver_cmd[@]}" > "${log_file}" 2>&1
  else
    OMP_NUM_THREADS=1 "${MPIEXEC}" -n "${np}" "${MPI_FLAG_ARR[@]}" \
      "${solver_cmd[@]}" > "${log_file}" 2>&1
  fi

  python3 - <<PY
import csv
import re
from pathlib import Path

log_path = Path("${log_file}")
csv_path = Path("${CSV}")
text = log_path.read_text()
runtime = re.search(r"Run time =\s*([0-9.eE+-]+)", text)
timing = re.search(
    r"Diag timing summary \\(max over MPI ranks\\): fillpatch=([0-9.eE+-]+) advance=([0-9.eE+-]+) reflux=([0-9.eE+-]+) tagging=([0-9.eE+-]+) total_runtime=([0-9.eE+-]+)",
    text,
)
fractions = re.search(
    r"Diag timing fractions \\(% of total runtime\\): fillpatch=([0-9.eE+-]+) advance=([0-9.eE+-]+) reflux=([0-9.eE+-]+) tagging=([0-9.eE+-]+)",
    text,
)
pf = re.search(
    r"Diag pressure-floor triggers \\(sum over MPI ranks\\): stage1=(\\d+) stage2=(\\d+) total=(\\d+)",
    text,
)
if not (runtime and timing and fractions and pf):
    raise RuntimeError(f"Could not parse diagnostic lines from {log_path}")

with csv_path.open("a", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(
        [
            "${case_id}",
            "${np}",
            runtime.group(1),
            timing.group(1),
            timing.group(2),
            timing.group(3),
            timing.group(4),
            fractions.group(1),
            fractions.group(2),
            fractions.group(3),
            fractions.group(4),
            pf.group(1),
            pf.group(2),
            pf.group(3),
        ]
    )
PY
}

run_case "task3_diag_np1" 1
run_case "task3_diag_np4" 4
run_case "task3_diag_np8" 8

python3 - <<PY
import csv
from pathlib import Path

csv_path = Path("${CSV}")
md_path = Path("${SUMMARY_MD}")
rows = list(csv.DictReader(csv_path.open()))

lines = []
lines.append("# Task 3 Diagnostic Timing Breakdown")
lines.append("")
lines.append("| Case | p | Runtime (s) | FillPatch (%) | Advance (%) | Reflux (%) | Tagging (%) | Pressure-floor triggers |")
lines.append("|---|---:|---:|---:|---:|---:|---:|---:|")
for row in rows:
    lines.append(
        f"| {row['case_id']} | {row['np']} | {float(row['walltime_sec']):.3f} | "
        f"{float(row['fillpatch_pct']):.2f} | {float(row['advance_pct']):.2f} | "
        f"{float(row['reflux_pct']):.2f} | {float(row['tagging_pct']):.2f} | "
        f"{row['pressure_floor_total']} |"
    )
md_path.write_text("\\n".join(lines) + "\\n")
PY

echo "Task 3 diagnostics written to ${RESULTS_DIR}"
