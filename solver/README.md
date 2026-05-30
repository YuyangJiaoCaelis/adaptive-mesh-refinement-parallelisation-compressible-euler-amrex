# Solver

This folder contains the final AMReX C++ solver source files and the minimal GNU Make build used by the experiment workflows.

## What is included

- `shared_solver/`: modified AMReX Euler solver source files
- `build/`: minimal AMReX build files
- `CHANGES_FROM_TUTORIAL.md`: provenance note describing which parts were adapted from the AMReX tutorial and which parts were written for this project

## Build Assumptions

- Linux is the recommended build target.
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

## Notes

- Experiment inputs, runner scripts, and analysis code live in `../experiments/`.
- AMReX plotting tools such as `fextract.gnu.ex` and `fcompare.gnu.ex` are expected under `AMREX_HOME/Tools/Plotfile/` by the experiment workflows.
