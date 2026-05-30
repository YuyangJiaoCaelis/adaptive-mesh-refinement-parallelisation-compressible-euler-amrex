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

RUN_SCOPE="${RUN_SCOPE:-all}" # smoke | full | all

N0=768
N1=1536
N2=3072

STOP_TIME=0.3
MPIEXEC="${MPIEXEC:-mpiexec}"
MPI_FLAGS="${MPI_FLAGS:---bind-to none --map-by :OVERSUBSCRIBE}"
read -r -a MPI_FLAG_ARR <<< "${MPI_FLAGS}"

RESULTS_DIR="${1:-${TASK_DIR}/results_task3_quadrant_$(date +%Y%m%d_%H%M)}"
mkdir -p "${RESULTS_DIR}"/{logs,raw,plots,tables}
MASTER_LOG="${RESULTS_DIR}/master.log"

TIME_UNIFORM_CSV="${RESULTS_DIR}/time_uniform.csv"
TIME_MPI_CSV="${RESULTS_DIR}/time_mpi.csv"
TIME_AMR_MATCH_CSV="${RESULTS_DIR}/time_amr_match.csv"
TIME_AMR_MPI_CSV="${RESULTS_DIR}/time_amr_mpi.csv"
SENS_CSV="${RESULTS_DIR}/sensitivity.csv"

SPEEDUP_CSV="${RESULTS_DIR}/speedup.csv"
EFF_CSV="${RESULTS_DIR}/efficiency.csv"
SUMMARY_MD="${RESULTS_DIR}/tables/task3_summary.md"

log_info() {
  echo "[INFO] $(date -Iseconds) $*" | tee -a "${MASTER_LOG}" >&2
}

log_warn() {
  echo "[WARN] $(date -Iseconds) $*" | tee -a "${MASTER_LOG}" >&2
}

append_csv_if_missing_header() {
  local csv="$1"
  local header="$2"
  if [[ ! -f "${csv}" ]]; then
    echo "${header}" > "${csv}"
  fi
}

append_csv_if_missing_header "${TIME_UNIFORM_CSV}" "N,cores,AMR_on,run_id,walltime_sec"
append_csv_if_missing_header "${TIME_MPI_CSV}" "N,cores,AMR_on,run_id,walltime_sec"
append_csv_if_missing_header "${TIME_AMR_MATCH_CSV}" "base_N,effective_N,max_level,cores,AMR_on,run_id,walltime_sec,level0_coverage_pct,level1_coverage_pct,level2_coverage_pct"
append_csv_if_missing_header "${TIME_AMR_MPI_CSV}" "base_N,effective_N,max_level,cores,AMR_on,run_id,walltime_sec,level0_coverage_pct,level1_coverage_pct,level2_coverage_pct"
append_csv_if_missing_header "${SENS_CSV}" "case_id,min_grid_size,n_error_buf,tagging_threshold,max_level,cores,run_id,walltime_sec,level0_coverage_pct,level1_coverage_pct,level2_coverage_pct"
append_csv_if_missing_header "${SPEEDUP_CSV}" "mode,N_or_effectiveN,cores,walltime_sec,speedup"
append_csv_if_missing_header "${EFF_CSV}" "mode,N_or_effectiveN,cores,walltime_sec,efficiency"

log_info "Task 3 runner started"
log_info "host=$(hostname)"
log_info "cwd=${TASK_DIR}"
log_info "results_dir=${RESULTS_DIR}"
log_info "run_scope=${RUN_SCOPE}"
log_info "mpi_flags=${MPI_FLAGS}"

to_cmd_string() {
  local out=""
  local arg
  for arg in "$@"; do
    out+=$(printf '%q ' "${arg}")
  done
  echo "${out}"
}

report_failure() {
  local failed_cmd="$1"
  local log_file="$2"
  local missing_output="$3"

  {
    echo "[ERROR] first failing command: ${failed_cmd}"
    echo "[ERROR] relevant log excerpt from: ${log_file}"
    if [[ -f "${log_file}" ]]; then
      tail -n 80 "${log_file}"
    else
      echo "[ERROR] log file not found"
    fi
    if [[ -n "${missing_output}" && ! -e "${missing_output}" ]]; then
      echo "[ERROR] missing output file: ${missing_output}"
    fi
  } | tee -a "${MASTER_LOG}" >&2

  exit 1
}

row_exists_uniform() {
  local n="$1"
  local cores="$2"
  local amr_on="$3"
  local run_id="$4"
  awk -F, -v n="${n}" -v c="${cores}" -v a="${amr_on}" -v r="${run_id}" \
    'NR>1 && $1==n && $2==c && $3==a && $4==r {found=1; exit} END{exit(found?0:1)}' "${TIME_UNIFORM_CSV}"
}

row_exists_mpi() {
  local n="$1"
  local cores="$2"
  local amr_on="$3"
  local run_id="$4"
  awk -F, -v n="${n}" -v c="${cores}" -v a="${amr_on}" -v r="${run_id}" \
    'NR>1 && $1==n && $2==c && $3==a && $4==r {found=1; exit} END{exit(found?0:1)}' "${TIME_MPI_CSV}"
}

row_exists_amr_match() {
  local base_n="$1"
  local eff_n="$2"
  local max_level="$3"
  local cores="$4"
  local amr_on="$5"
  local run_id="$6"
  awk -F, -v b="${base_n}" -v e="${eff_n}" -v l="${max_level}" -v c="${cores}" -v a="${amr_on}" -v r="${run_id}" \
    'NR>1 && $1==b && $2==e && $3==l && $4==c && $5==a && $6==r {found=1; exit} END{exit(found?0:1)}' "${TIME_AMR_MATCH_CSV}"
}

row_exists_amr_mpi() {
  local base_n="$1"
  local eff_n="$2"
  local max_level="$3"
  local cores="$4"
  local amr_on="$5"
  local run_id="$6"
  awk -F, -v b="${base_n}" -v e="${eff_n}" -v l="${max_level}" -v c="${cores}" -v a="${amr_on}" -v r="${run_id}" \
    'NR>1 && $1==b && $2==e && $3==l && $4==c && $5==a && $6==r {found=1; exit} END{exit(found?0:1)}' "${TIME_AMR_MPI_CSV}"
}

row_exists_sens() {
  local case_id="$1"
  awk -F, -v id="${case_id}" 'NR>1 && $1==id {found=1; exit} END{exit(found?0:1)}' "${SENS_CSV}"
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

lookup_uniform_runtime() {
  local n="$1"
  local run_id="$2"
  awk -F, -v n="${n}" -v r="${run_id}" 'NR>1 && $1==n && $2==1 && $3==0 && $4==r {print $5; exit}' "${TIME_UNIFORM_CSV}"
}

lookup_amr_match_row() {
  local base_n="$1"
  local eff_n="$2"
  local max_level="$3"
  local run_id="$4"
  awk -F, -v b="${base_n}" -v e="${eff_n}" -v l="${max_level}" -v r="${run_id}" \
    'NR>1 && $1==b && $2==e && $3==l && $4==1 && $5==1 && $6==r {print $7","$8","$9","$10; exit}' "${TIME_AMR_MATCH_CSV}"
}

run_solver_case() {
  local case_id="$1"
  local np="$2"
  local base_n="$3"
  local max_level="$4"
  local regrid_int="$5"
  local nbuf="$6"
  local blocking_factor="$7"
  local max_grid_size="$8"
  local phigrad="$9"
  local do_reflux=0
  if [[ "${max_level}" -gt 0 ]]; then
    do_reflux=1
  fi

  local log_file="${RESULTS_DIR}/logs/${case_id}.log"
  local plot_prefix="${RESULTS_DIR}/raw/plt_${case_id}"

  rm -rf "${plot_prefix}"*

  local -a cmd=(
    "${EXE}" inputs
    "prob.ic=3"
    "prob.x0=0.5"
    "prob.y0=0.5"
    "stop_time=${STOP_TIME}"
    "adv.first_order=0"
    "adv.force_1d=0"
    "adv.host_update=0"
    "adv.debug_state=0"
    "adv.v=0"
    "amr.v=0"
    "adv.do_reflux=${do_reflux}"
    "geometry.is_periodic=0 0"
    "amr.max_level=${max_level}"
    "amr.n_cell=${base_n} ${base_n}"
    "amr.ref_ratio=2 2 2 2"
    "amr.regrid_int=${regrid_int}"
    "amr.n_error_buf=${nbuf}"
    "amr.blocking_factor=${blocking_factor}"
    "amr.max_grid_size=${max_grid_size}"
    "amr.grid_eff=0.75"
    "amr.plot_int=1000000"
    "amr.plot_file=${plot_prefix}"
    "tagging.max_phierr_lev=-1"
    "tagging.max_phigrad_lev=${max_level}"
    "tagging.phigrad=${phigrad}"
  )

  log_info "RUN case=${case_id} np=${np} base_n=${base_n} max_level=${max_level}"

  if [[ "${np}" -eq 1 ]]; then
    local cmd_str
    cmd_str="$(to_cmd_string "${cmd[@]}")"
    echo "[CMD] ${cmd_str}" >> "${MASTER_LOG}"
    if ! OMP_NUM_THREADS=1 "${cmd[@]}" > "${log_file}" 2>&1; then
      report_failure "${cmd_str}" "${log_file}" "${plot_prefix}"
    fi
  else
    local -a mpi_cmd=("${MPIEXEC}" -n "${np}" "${MPI_FLAG_ARR[@]}" "${cmd[@]}")
    local mpi_cmd_str
    mpi_cmd_str="$(to_cmd_string "${mpi_cmd[@]}")"
    echo "[CMD] ${mpi_cmd_str}" >> "${MASTER_LOG}"
    if ! "${mpi_cmd[@]}" > "${log_file}" 2>&1; then
      report_failure "${mpi_cmd_str}" "${log_file}" "${plot_prefix}"
    fi
  fi

  local rt
  rt="$(parse_runtime "${log_file}")"
  if [[ -z "${rt}" ]]; then
    report_failure "parse runtime from ${log_file}" "${log_file}" "${log_file}"
  fi

  local cov
  cov="$(parse_coverages "${log_file}")"
  local l0 l1 l2
  IFS=',' read -r l0 l1 l2 <<< "${cov}"

  log_info "DONE case=${case_id} runtime_s=${rt} coverage=(${l0},${l1},${l2})"
  echo "${rt},${l0},${l1},${l2}"
}

run_uniform_serial_case() {
  local n="$1"
  local run_id="$2"

  if row_exists_uniform "${n}" 1 0 "${run_id}"; then
    log_info "SKIP existing uniform serial N=${n} run_id=${run_id}"
    return
  fi

  local case_id="uniform_n${n}_np1_${run_id}"
  local out
  out="$(run_solver_case "${case_id}" 1 "${n}" 0 2 2 16 64 "0.08 0.04")"
  local rt
  IFS=',' read -r rt _ <<< "${out}"
  echo "${n},1,0,${run_id},${rt}" >> "${TIME_UNIFORM_CSV}"
}

copy_uniform_to_mpi_case() {
  local n="$1"
  local run_id="$2"
  local rt
  rt="$(lookup_uniform_runtime "${n}" "${run_id}")"
  if [[ -z "${rt}" ]]; then
    report_failure "lookup uniform runtime for N=${n}, run_id=${run_id}" "${MASTER_LOG}" "${TIME_UNIFORM_CSV}"
  fi
  echo "${n},1,0,${run_id},${rt}" >> "${TIME_MPI_CSV}"
  log_info "COPIED p=1 uniform timing into time_mpi.csv for N=${n}, run_id=${run_id}, rt=${rt}"
}

run_uniform_mpi_case() {
  local n="$1"
  local np="$2"
  local run_id="$3"

  if row_exists_mpi "${n}" "${np}" 0 "${run_id}"; then
    log_info "SKIP existing uniform MPI N=${n} np=${np} run_id=${run_id}"
    return
  fi

  if [[ "${np}" -eq 1 ]]; then
    # Reuse part-1 serial timings for p=1 to avoid duplicate expensive runs.
    copy_uniform_to_mpi_case "${n}" "${run_id}"
    return
  fi

  local case_id="uniform_n${n}_np${np}_${run_id}"
  local out
  out="$(run_solver_case "${case_id}" "${np}" "${n}" 0 2 2 16 64 "0.08 0.04")"
  local rt
  IFS=',' read -r rt _ <<< "${out}"
  echo "${n},${np},0,${run_id},${rt}" >> "${TIME_MPI_CSV}"
}

run_amr_match_case() {
  local max_level="$1"
  local effective_n="$2"
  local run_id="$3"
  local phigrad="$4"

  if row_exists_amr_match "${N0}" "${effective_n}" "${max_level}" 1 1 "${run_id}"; then
    log_info "SKIP existing AMR match effective_N=${effective_n} max_level=${max_level} run_id=${run_id}"
    return
  fi

  local case_id="amr_match_eff${effective_n}_np1_${run_id}"
  local out
  out="$(run_solver_case "${case_id}" 1 "${N0}" "${max_level}" 2 2 16 64 "${phigrad}")"
  local rt l0 l1 l2
  IFS=',' read -r rt l0 l1 l2 <<< "${out}"
  echo "${N0},${effective_n},${max_level},1,1,${run_id},${rt},${l0},${l1},${l2}" >> "${TIME_AMR_MATCH_CSV}"
}

copy_amr_match_to_amr_mpi_p1() {
  local run_id="$1"
  local ref_row
  ref_row="$(lookup_amr_match_row "${N0}" "${N2}" 2 "${run_id}")"
  if [[ -z "${ref_row}" ]]; then
    report_failure "lookup AMR match row for highest effective resolution" "${MASTER_LOG}" "${TIME_AMR_MATCH_CSV}"
  fi
  echo "${N0},${N2},2,1,1,${run_id},${ref_row}" >> "${TIME_AMR_MPI_CSV}"
  log_info "COPIED AMR match p=1 timing into time_amr_mpi.csv for effective_N=${N2}, run_id=${run_id}"
}

run_amr_mpi_high_case() {
  local np="$1"
  local run_id="$2"

  if row_exists_amr_mpi "${N0}" "${N2}" 2 "${np}" 1 "${run_id}"; then
    log_info "SKIP existing AMR MPI high np=${np} run_id=${run_id}"
    return
  fi

  if [[ "${np}" -eq 1 ]]; then
    # Reuse part-3 p=1 AMR-high run.
    copy_amr_match_to_amr_mpi_p1 "${run_id}"
    return
  fi

  local case_id="amr_high_eff${N2}_np${np}_${run_id}"
  local out
  out="$(run_solver_case "${case_id}" "${np}" "${N0}" 2 2 2 16 64 "0.08 0.04")"
  local rt l0 l1 l2
  IFS=',' read -r rt l0 l1 l2 <<< "${out}"
  echo "${N0},${N2},2,${np},1,${run_id},${rt},${l0},${l1},${l2}" >> "${TIME_AMR_MPI_CSV}"
}

run_sensitivity_case() {
  local case_id="$1"
  local min_grid_size="$2"
  local nbuf="$3"
  local threshold="$4"
  local run_id="$5"

  if row_exists_sens "${case_id}"; then
    log_info "SKIP existing sensitivity case=${case_id}"
    return
  fi

  local max_grid_size=64
  if [[ "${min_grid_size}" -gt "${max_grid_size}" ]]; then
    max_grid_size="${min_grid_size}"
  fi

  local full_case_id="sens_${case_id}_${run_id}"
  local out
  out="$(run_solver_case "${full_case_id}" 1 "${N0}" 1 2 "${nbuf}" "${min_grid_size}" "${max_grid_size}" "${threshold}")"
  local rt l0 l1 l2
  IFS=',' read -r rt l0 l1 l2 <<< "${out}"
  echo "${case_id},${min_grid_size},${nbuf},${threshold},1,1,${run_id},${rt},${l0},${l1},${l2}" >> "${SENS_CSV}"
}

run_smoke() {
  log_info "=== Smoke phase: Part (1) N0/N1 once ==="
  run_uniform_serial_case "${N0}" r1
  run_uniform_serial_case "${N1}" r1
}

run_full() {
  log_info "=== Full phase: Task 3 full matrix ==="

  # Part (1): uniform no-AMR single core.
  run_uniform_serial_case "${N0}" r1
  run_uniform_serial_case "${N0}" r2
  run_uniform_serial_case "${N0}" r3
  run_uniform_serial_case "${N1}" r1
  run_uniform_serial_case "${N1}" r2
  run_uniform_serial_case "${N2}" r1

  # Part (2): uniform no-AMR MPI scaling.
  for np in 1 2 4 8 16; do
    for run_id in r1 r2; do
      run_uniform_mpi_case "${N0}" "${np}" "${run_id}"
    done
  done
  for np in 1 2 4 8 16; do
    run_uniform_mpi_case "${N1}" "${np}" r1
  done
  for np in 1 4 8 16; do
    run_uniform_mpi_case "${N2}" "${np}" r1
  done

  # Part (3): AMR effective-resolution matching.
  run_amr_match_case 1 "${N1}" r1 "0.08"
  run_amr_match_case 2 "${N2}" r1 "0.08 0.04"

  # Part (4): AMR+MPI high effective resolution.
  for np in 1 2 4 8 16; do
    run_amr_mpi_high_case "${np}" r1
  done

  # Part (5): sensitivity one-at-a-time around baseline (32, 2, 0.1).
  run_sensitivity_case "baseline" 32 2 0.1 r1
  run_sensitivity_case "min_grid_16" 16 2 0.1 r1
  run_sensitivity_case "min_grid_64" 64 2 0.1 r1
  run_sensitivity_case "n_error_buf_1" 32 1 0.1 r1
  run_sensitivity_case "n_error_buf_4" 32 4 0.1 r1
  run_sensitivity_case "threshold_005" 32 2 0.05 r1
  run_sensitivity_case "threshold_020" 32 2 0.2 r1
}

log_info "Phases requested: ${RUN_SCOPE}"

case "${RUN_SCOPE}" in
  smoke)
    run_smoke
    ;;
  full)
    run_full
    ;;
  all)
    run_smoke
    run_full
    ;;
  *)
    echo "ERROR: unknown RUN_SCOPE=${RUN_SCOPE} (use smoke|full|all)" | tee -a "${MASTER_LOG}" >&2
    exit 1
    ;;
esac

ANALYSIS_CMD=(
  python3 "${TASK_DIR}/analysis/summarize_task3_quadrant.py"
  --results-dir "${RESULTS_DIR}"
)
log_info "Running post-processing summary"
echo "[CMD] $(to_cmd_string "${ANALYSIS_CMD[@]}")" >> "${MASTER_LOG}"
if ! "${ANALYSIS_CMD[@]}" >> "${MASTER_LOG}" 2>&1; then
  report_failure "$(to_cmd_string "${ANALYSIS_CMD[@]}")" "${MASTER_LOG}" "${SUMMARY_MD}"
fi

log_info "Task 3 run complete. Summary: ${SUMMARY_MD}"
