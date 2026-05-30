# Experiments

This folder contains the reproducible workflows behind the report results.

## Layout

- `riemann_1d/`: Toro shock-tube validation, AMR comparisons, and smooth convergence checks.
- `orientation_2d/`: two-dimensional orientation consistency and diagonal split tests.
- `quadrant_scaling/`: Lax-Liu-style quadrant benchmark for AMR localisation, MPI scaling, and full-field accuracy diagnostics.

Each experiment has `inputs/`, `scripts/`, and `analysis/` subdirectories. Generated result directories are ignored by git; the frozen report outputs live under `../results/`.
