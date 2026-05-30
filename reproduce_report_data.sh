#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")" && pwd)"

make -C "${ROOT_DIR}" build

"${ROOT_DIR}/experiments/riemann_1d/scripts/run_toro1d_uniform_amr.sh"
"${ROOT_DIR}/experiments/riemann_1d/scripts/run_smooth_entropy_convergence.sh"

"${ROOT_DIR}/experiments/orientation_2d/scripts/run_task2_2d_fullmatrix_report.sh"

RUN_SCOPE="${RUN_SCOPE:-all}" \
  "${ROOT_DIR}/experiments/quadrant_scaling/scripts/run_task3_quadrant_matrix.sh"
"${ROOT_DIR}/experiments/quadrant_scaling/scripts/run_task3_highres_refresh.sh"
"${ROOT_DIR}/experiments/quadrant_scaling/scripts/run_task3_accuracy_refresh.sh"
"${ROOT_DIR}/experiments/quadrant_scaling/scripts/run_task3_diagnostics.sh"
