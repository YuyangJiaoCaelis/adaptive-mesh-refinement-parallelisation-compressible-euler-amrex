# Changes From Tutorial

Modified content in this repository is attributable to this public project version.

This solver package is based on the AMReX advection tutorial used as the recommended starting point in the assignment brief.

## Tutorial base

The main AMR framework files were adapted from:

- `amrex-tutorials/ExampleCodes/Amr/Advection_AmrLevel/Source/AmrLevelAdv.cpp`
- `amrex-tutorials/ExampleCodes/Amr/Advection_AmrLevel/Source/AmrLevelAdv.H`
- `amrex-tutorials/ExampleCodes/Amr/Advection_AmrLevel/Source/Adv.cpp`
- `amrex-tutorials/ExampleCodes/Amr/Advection_AmrLevel/Exec/UniformVelocity/Adv_prob.cpp`
- `amrex-tutorials/ExampleCodes/Amr/Advection_AmrLevel/Exec/UniformVelocity/Prob_Parm.H`

This repository keeps only the files needed to reproduce the data quoted in the report.

## Main solver changes

### `shared_solver/Adv.cpp`
Adapted from the tutorial advection update into a 2D compressible Euler solver. The main changes are:

- replaced scalar advection with conservative Euler variables `(rho, momx, momy, E)`
- added primitive-variable reconstruction and conversion between conservative and primitive states
- added HLL fluxes in both spatial directions
- added SSPRK2 stage updates
- added positivity floors for density and pressure
- added the stage-wise pressure-floor energy reset used in the report
- added optional timing and pressure-floor diagnostics used for the analysis

### `shared_solver/AmrLevelAdv.cpp`
Adapted from the tutorial AMR level driver. The main changes are:

- changed the state definition from a single scalar to Euler conserved variables
- changed boundary-condition handling from periodic-only tutorial defaults to the periodic/outflow combinations used in the report
- added reflux and AMR diagnostic handling used in the timing study
- added pressure-floor diagnostic counters
- added the task-specific tagging logic and timing bookkeeping used in the report

### `shared_solver/AmrLevelAdv.H`
Adapted from the tutorial header. The main changes are:

- updated the state layout for Euler variables instead of a scalar field
- added boundary-condition storage
- added diagnostic members for timing and pressure-floor counts
- trimmed inherited tutorial comments that were not needed for this project

### `shared_solver/Prob.cpp`
This file is project-specific rather than a direct tutorial copy. It defines the problem initial data used in the report:

- Toro shock-tube tests
- axis-aligned and diagonal 2D extensions
- Lax-Liu quadrant problem
- smooth entropy-wave test cases

### `shared_solver/Prob.H` and `shared_solver/Prob_Parm.H`
These files were simplified for the Euler setup. The tutorial advection-velocity parameter structure is no longer used.

### `shared_solver/Adv_prob.cpp`
The tutorial `probin` initialisation logic for constant advection velocity was removed because the Euler cases in this project do not use that mechanism.

## Build files

### `build/GNUmakefile`, `build/Make.Adv`, `build/Make.package`
These are trimmed build files needed to compile the included solver sources against AMReX.

## Task-specific scripts and analysis

The tutorial does not provide the task runners or post-processing used for this report. The following directories were written for this project to reproduce the quoted report data:

- `experiments/riemann_1d/`
- `experiments/orientation_2d/`
- `experiments/quadrant_scaling/`

These contain:

- input files used for the quoted runs
- shell scripts that run the specific report matrices
- Python analysis scripts used to generate the quoted error tables, timing summaries, sensitivity summaries, and full-field checks

## Scope of this pack

This is not a full development tree. It contains only:

- the modified solver files needed for the quoted results
- the input files used in the quoted runs
- the run scripts and analysis scripts needed to regenerate the quoted report data

The original AMReX tutorial and the full working project tree are not included in this submission folder.
