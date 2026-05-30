#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
TASK_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
CODE_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
cd "${TASK_DIR}"

EXE="${CODE_ROOT}/build/main2d.gnu.MPI.ex"
RESULTS_DIR="${1:-${TASK_DIR}/results_pressure_floor_representatives_$(date +%Y%m%d_%H%M%S)}"
mkdir -p "${RESULTS_DIR}"/{logs,raw,checks}

CSV="${RESULTS_DIR}/checks/pressure_floor_representatives.csv"
MD="${RESULTS_DIR}/checks/pressure_floor_representatives.md"
echo "case_id,walltime_sec,pressure_floor_stage1,pressure_floor_stage2,pressure_floor_total" > "${CSV}"

run_case() {
  local case_id="$1"
  shift
  local log_file="${RESULTS_DIR}/logs/${case_id}.log"
  local plot_prefix="${RESULTS_DIR}/raw/plt_${case_id}"
  rm -rf "${plot_prefix}"*

  OMP_NUM_THREADS=1 "${EXE}" inputs \
    "$@" \
    "adv.diag_pressure_floor=1" \
    "adv.v=0" \
    "amr.v=0" \
    "amr.plot_int=1000000" \
    "amr.plot_file=${plot_prefix}" \
    > "${log_file}" 2>&1

  python3 - <<PY
import csv
import re
from pathlib import Path

log_path = Path("${log_file}")
csv_path = Path("${CSV}")
text = log_path.read_text()
runtime = re.search(r"Run time =\s*([0-9.eE+-]+)", text)
pf = re.search(
    r"Diag pressure-floor triggers \\(sum over MPI ranks\\): stage1=(\\d+) stage2=(\\d+) total=(\\d+)",
    text,
)
if not (runtime and pf):
    raise RuntimeError(f"Could not parse diagnostics from {log_path}")
with csv_path.open("a", newline="") as f:
    csv.writer(f).writerow(
        ["${case_id}", runtime.group(1), pf.group(1), pf.group(2), pf.group(3)]
    )
PY
}

run_case \
  "task1_test5_amr" \
  "prob.test=5" \
  "prob.ic=0" \
  "prob.print_ic=0" \
  "stop_time=0.035" \
  "prob.x0=0.5" \
  "adv.first_order=0" \
  "adv.force_1d=1" \
  "adv.host_update=0" \
  "adv.debug_state=0" \
  "adv.do_reflux=1" \
  "amr.max_level=2" \
  "amr.n_cell=100 4" \
  "amr.ref_ratio=2 2 2 2" \
  "amr.regrid_int=2" \
  "amr.n_error_buf=2" \
  "amr.blocking_factor=4" \
  "amr.max_grid_size=64" \
  "amr.grid_eff=0.7" \
  "geometry.is_periodic=0 1" \
  "tagging.max_phierr_lev=-1" \
  "tagging.max_phigrad_lev=2" \
  "tagging.phigrad=0.04 0.02"

run_case \
  "task2_test4_diag_amr" \
  "prob.test=4" \
  "prob.ic=2" \
  "prob.print_ic=0" \
  "prob.x0=0.5" \
  "prob.y0=0.5" \
  "prob.diag_nx=1.0" \
  "prob.diag_ny=1.0" \
  "stop_time=0.035" \
  "adv.first_order=0" \
  "adv.force_1d=0" \
  "adv.host_update=0" \
  "adv.debug_state=0" \
  "adv.do_reflux=1" \
  "geometry.is_periodic=0 0" \
  "amr.max_level=1" \
  "amr.n_cell=200 200" \
  "amr.ref_ratio=2 2 2 2" \
  "amr.regrid_int=2" \
  "amr.n_error_buf=2" \
  "amr.blocking_factor=8" \
  "amr.max_grid_size=64" \
  "amr.grid_eff=0.7" \
  "tagging.max_phierr_lev=-1" \
  "tagging.max_phigrad_lev=1" \
  "tagging.phigrad=0.06 0.03"

python3 - <<PY
import csv
from pathlib import Path

csv_path = Path("${CSV}")
md_path = Path("${MD}")
rows = list(csv.DictReader(csv_path.open()))
lines = []
lines.append("# Representative Pressure-Floor Diagnostics")
lines.append("")
lines.append("| Case | Runtime (s) | Stage-1 triggers | Stage-2 triggers | Total |")
lines.append("|---|---:|---:|---:|---:|")
for row in rows:
    lines.append(
        f"| {row['case_id']} | {float(row['walltime_sec']):.3f} | "
        f"{row['pressure_floor_stage1']} | {row['pressure_floor_stage2']} | {row['pressure_floor_total']} |"
    )
md_path.write_text("\\n".join(lines) + "\\n")
PY

echo "Pressure-floor diagnostics written to ${RESULTS_DIR}"
