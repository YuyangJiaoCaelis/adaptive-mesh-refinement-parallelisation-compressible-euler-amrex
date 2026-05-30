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

if [[ $# -ge 1 ]]; then
  RESULTS_DIR="$1"
else
  TS="$(date +%Y%m%d_%H%M)"
  RESULTS_DIR="${TASK_DIR}/results_task2_2d_fullmatrix_${TS}"
fi

mkdir -p "${RESULTS_DIR}"/{logs,plotfiles,figures/main_text,figures/appendix,checks}
printf "%s\n" "${RESULTS_DIR}" > "${TASK_DIR}/.last_task2_fullmatrix_result_dir"

# Numerical settings aligned with RUN_MATRIX.md task-2 recommendations.
N_UNIFORM=400
N_AMR_BASE=200
MAX_LEVEL_AMR=1

declare -A STOP_TIME
STOP_TIME[1]=0.25
STOP_TIME[2]=0.15
STOP_TIME[3]=0.012
STOP_TIME[4]=0.035
STOP_TIME[5]=0.035

declare -A SPLIT_NAME
SPLIT_NAME[0]=x
SPLIT_NAME[1]=y
SPLIT_NAME[2]=diag

latest_plotfile() {
  local prefix="$1"
  python3 - <<'PY' "$prefix"
import pathlib, re, sys
prefix = sys.argv[1]
p = pathlib.Path(prefix)
parent = p.parent
stem = p.name
pat = re.compile(rf"^{re.escape(stem)}\d+$")
cands = sorted([x for x in parent.glob(f"{stem}*") if x.is_dir() and pat.match(x.name)])
print(cands[-1] if cands else "")
PY
}

run_case() {
  local test_id="$1"
  local ic="$2"
  local mode="$3" # uniform|amr

  local stop_time="${STOP_TIME[${test_id}]}"
  local split="${SPLIT_NAME[${ic}]}"
  local max_level
  local ncell
  local do_reflux
  local tag_lev

  if [[ "${mode}" == "uniform" ]]; then
    max_level=0
    ncell="${N_UNIFORM} ${N_UNIFORM}"
    do_reflux=0
    tag_lev=0
  else
    max_level="${MAX_LEVEL_AMR}"
    ncell="${N_AMR_BASE} ${N_AMR_BASE}"
    do_reflux=1
    tag_lev="${MAX_LEVEL_AMR}"
  fi

  local prefix="${RESULTS_DIR}/plotfiles/plt_t${test_id}_ic${ic}_${mode}"
  local runlog="${RESULTS_DIR}/logs/run_t${test_id}_ic${ic}_${mode}.log"

  echo "[RUN] test=${test_id} ic=${ic} (${split}) mode=${mode}"

  "${EXE}" inputs \
    "prob.test=${test_id}" \
    "prob.ic=${ic}" \
    "prob.x0=0.5" \
    "prob.y0=0.5" \
    "prob.diag_nx=1.0" \
    "prob.diag_ny=1.0" \
    "stop_time=${stop_time}" \
    "adv.first_order=0" \
    "adv.force_1d=0" \
    "adv.host_update=0" \
    "adv.debug_state=0" \
    "adv.do_reflux=${do_reflux}" \
    "geometry.is_periodic=0 0" \
    "amr.max_level=${max_level}" \
    "amr.n_cell=${ncell}" \
    "amr.ref_ratio=2 2 2 2" \
    "amr.regrid_int=2" \
    "amr.n_error_buf=2" \
    "amr.blocking_factor=8" \
    "amr.max_grid_size=64" \
    "amr.grid_eff=0.7" \
    "amr.plot_int=1000000" \
    "amr.plot_file=${prefix}" \
    "tagging.max_phierr_lev=-1" \
    "tagging.max_phigrad_lev=${tag_lev}" \
    "tagging.phigrad=0.06 0.03" \
    > "${runlog}" 2>&1

  local plt
  plt="$(latest_plotfile "${prefix}")"
  if [[ -z "${plt}" ]]; then
    echo "ERROR: no plotfile found for test=${test_id} ic=${ic} mode=${mode}" >&2
    echo "Log excerpt:" >&2
    tail -n 40 "${runlog}" >&2 || true
    exit 2
  fi

  printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
    "${test_id}" "${ic}" "${split}" "${mode}" "${plt}" "${runlog}" \
    >> "${RESULTS_DIR}/logs/run_index.tsv"

  echo "[DONE] test=${test_id} ic=${ic} mode=${mode} -> $(basename "${plt}")"
}

{
  printf "test\tic\tsplit\tmode\tplotfile\tlog\n"
} > "${RESULTS_DIR}/logs/run_index.tsv"

for test_id in 1 2 3 4 5; do
  for ic in 0 1 2; do
    run_case "${test_id}" "${ic}" uniform
    run_case "${test_id}" "${ic}" amr
  done
done

cp -f "${TASK_DIR}/inputs" "${RESULTS_DIR}/inputs.base"
cat > "${RESULTS_DIR}/checks/resolution_note.txt" <<'EOF'
Task-2 2D matched-resolution setup:
- Uniform: amr.n_cell = 400 400, amr.max_level = 0
- AMR:     amr.n_cell = 200 200, amr.max_level = 1, amr.ref_ratio = 2
- Effective finest spacing is matched at approximately 1/400 in each direction.
- amr.blocking_factor = 8 divides both 400 and 200 exactly.
EOF

python3 "${TASK_DIR}/analysis/analyze_task2_2d_fullmatrix.py" --results-dir "${RESULTS_DIR}" \
  > "${RESULTS_DIR}/logs/analysis.log" 2>&1

echo "All runs and analysis completed."
echo "RESULTS_DIR=${RESULTS_DIR}"
