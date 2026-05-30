#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
TASK_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
CODE_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
cd "${TASK_DIR}"

EXE="${CODE_ROOT}/build/main2d.gnu.MPI.ex"
if [[ ! -x "${EXE}" ]]; then
  echo "ERROR: executable not found: ${EXE}" >&2
  echo "Build first with: make -j8" >&2
  exit 1
fi

if [[ -z "${AMREX_HOME:-}" ]]; then
  echo "ERROR: set AMREX_HOME so fextract.gnu.ex can be found." >&2
  exit 1
fi
FEXTRACT="${AMREX_HOME}/Tools/Plotfile/fextract.gnu.ex"
if [[ ! -x "${FEXTRACT}" ]]; then
  echo "ERROR: fextract not found: ${FEXTRACT}" >&2
  exit 1
fi

RESULTS_DIR="${1:-${TASK_DIR}/results_toro1d_uniform_vs_amr}"
mkdir -p "${RESULTS_DIR}"

declare -A STOP_TIME
STOP_TIME[1]=0.25
STOP_TIME[2]=0.15
STOP_TIME[3]=0.012
STOP_TIME[4]=0.035
STOP_TIME[5]=0.035

write_inputs_file() {
  local test_id="$1"
  local stop_time="$2"
  local out="${RESULTS_DIR}/inputs.t${test_id}"
  cat > "${out}" <<EOF
stop_time = ${stop_time}
prob.test = ${test_id}
prob.gamma = 1.4
prob.x0 = 0.5
EOF
  echo "${out}"
}

latest_plotfile() {
  local prefix="$1"
  local latest
  latest="$(ls -d "${prefix}"* 2>/dev/null | grep -E "${prefix}[0-9]+$" | sort | tail -n 1 || true)"
  echo "${latest}"
}

run_case() {
  local mode="$1"      # uniform or amr
  local test_id="$2"
  local eff_n="$3"     # 100/200/400

  local stop_time="${STOP_TIME[${test_id}]}"
  local max_level
  local base_n
  local phigrad
  local do_reflux

  if [[ "${mode}" == "uniform" ]]; then
    max_level=0
    base_n="${eff_n}"
    phigrad="0.04 0.02"
    do_reflux=0
  else
    case "${eff_n}" in
      100) max_level=0 ;;
      200) max_level=1 ;;
      400) max_level=2 ;;
      *) echo "ERROR: unsupported effective N=${eff_n} for AMR" >&2; exit 2 ;;
    esac
    base_n=100
    phigrad="0.04 0.02"
    do_reflux=1
  fi

  local prefix="${RESULTS_DIR}/plt_${mode}_t${test_id}_n${eff_n}"
  local runlog="${RESULTS_DIR}/${mode}_t${test_id}_n${eff_n}.log"
  local pplog="${RESULTS_DIR}/${mode}_t${test_id}_n${eff_n}_pp.log"
  local input_file
  input_file="$(write_inputs_file "${test_id}" "${stop_time}")"

  rm -rf "${prefix}"*

  "${EXE}" inputs \
    "prob.test=${test_id}" \
    "prob.ic=0" \
    "prob.print_ic=1" \
    "stop_time=${stop_time}" \
    "prob.x0=0.5" \
    "adv.first_order=0" \
    "adv.force_1d=1" \
    "adv.host_update=0" \
    "adv.debug_state=0" \
    "adv.do_reflux=${do_reflux}" \
    "amr.max_level=${max_level}" \
    "amr.n_cell=${base_n} 4" \
    "amr.ref_ratio=2 2 2 2" \
    "amr.regrid_int=2" \
    "amr.n_error_buf=2" \
    "amr.blocking_factor=4" \
    "amr.max_grid_size=64" \
    "amr.grid_eff=0.7" \
    "amr.plot_int=1000000" \
    "amr.plot_file=${prefix}" \
    "geometry.is_periodic=0 1" \
    "tagging.max_phierr_lev=-1" \
    "tagging.max_phigrad_lev=${max_level}" \
    "tagging.phigrad=${phigrad}" \
    > "${runlog}" 2>&1

  local plt
  plt="$(latest_plotfile "${prefix}")"
  if [[ -z "${plt}" ]]; then
    echo "ERROR: no plotfile for ${mode} t${test_id} n${eff_n}" >&2
    exit 3
  fi

  local derived="${RESULTS_DIR}/${mode}_t${test_id}_n${eff_n}_derived.csv"
  local exact="${RESULTS_DIR}/${mode}_t${test_id}_n${eff_n}_exact.csv"

  if [[ "${mode}" == "uniform" ]]; then
    python3 postprocess_1d.py \
      --plotfile "${plt}" \
      --out "${derived}" \
      --inputs "${input_file}" \
      --exact "${exact}" \
      > "${pplog}" 2>&1
  else
    local slice="${RESULTS_DIR}/${mode}_t${test_id}_n${eff_n}.slice.dat"
    "${FEXTRACT}" -s "${slice}" -v "rho momx momy E" "${plt}" > /dev/null 2>&1
    python3 postprocess_1d.py \
      --slice "${slice}" \
      --out "${derived}" \
      --inputs "${input_file}" \
      --exact "${exact}" \
      --time "${stop_time}" \
      > "${pplog}" 2>&1
  fi

  echo "DONE ${mode} test=${test_id} N=${eff_n} -> ${plt}"
}

for mode in uniform amr; do
  for t in 1 2 3 4 5; do
    for n in 100 200 400; do
      run_case "${mode}" "${t}" "${n}"
    done
  done
done

python3 "${SCRIPT_DIR}/toro1d_error_table.py" \
  --results-dir "${RESULTS_DIR}" \
  --out-csv "toro1d_errors.csv" \
  --out-md "toro1d_errors.md"

echo "All task-1 runs complete. Results in ${RESULTS_DIR}"
