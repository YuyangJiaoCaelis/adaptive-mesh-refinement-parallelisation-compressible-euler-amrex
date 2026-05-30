# Code Submission README

This code folder contains the modified source files, run scripts, and analysis scripts used to reproduce the data quoted in the report for the project report.

## What is included

- `shared_solver/`: modified AMReX Euler solver source files
- `build/`: minimal AMReX build files
- `task1/`, `task2/`, `task3/`: task-specific input files, run scripts, and post-processing scripts used for the report
- `CHANGES_FROM_TUTORIAL.md`: provenance note describing which parts were adapted from the AMReX tutorial and which parts were written for this project
- `README_REPORT_SETTINGS.md`: note explaining that the task scripts override the base input templates where needed

## Build assumptions

- The assessors may assume Linux.
- AMReX is required and should be available through the `AMREX_HOME` environment variable.
- The build uses GNU Make and the trimmed files in `build/`.

Example build sequence:

```bash
cd build
make -j8 AMREX_HOME=/path/to/amrex
```

This produces the solver executable inside the build directory:

```bash
build/main2d.gnu.MPI.ex
```

## Running the quoted report workflows

The report settings are reproduced by the task runner scripts:

- `task1/scripts/run_toro1d_uniform_amr.sh`
- `task1/scripts/run_smooth_entropy_convergence.sh`
- `task1/scripts/run_smooth_amr_convergence.sh`
- `task1/scripts/run_smooth_multilevel_periodic_convergence.sh`
- `task2/scripts/run_task2_2d_fullmatrix_report.sh`
- `task3/scripts/run_task3_quadrant_matrix.sh`
- `task3/scripts/run_task3_highres_refresh.sh`
- `task3/scripts/run_task3_lowerres_refresh.sh`
- `task3/scripts/run_task3_accuracy_refresh.sh`
- `task3/scripts/run_task3_diagnostics.sh`
- `task3/scripts/run_pressure_floor_representatives.sh`

Each script writes a self-contained results directory and then calls the matching analysis scripts where required.

## Notes for assessors

- The task scripts set the report settings directly.
- The base templates in `task1/inputs`, `task2/inputs`, and `task3/inputs` are included for reference; the task scripts override them where needed.
- AMReX plotting tools such as `fextract.gnu.ex` and `fcompare.gnu.ex` are expected under `AMREX_HOME/Tools/Plotfile/`.
- Only the files needed to reproduce the quoted report data are included in this pack.
