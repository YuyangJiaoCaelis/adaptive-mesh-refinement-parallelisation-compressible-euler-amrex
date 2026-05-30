# Numerical Method Evidence (Task 1)

## Solver class and order configuration

- The Task-1 run script enforces second-order mode by setting `adv.first_order=0` in every run command:
  - `scripts/run_toro1d_uniform_amr.sh` line 97.
- The same script sets `prob.gamma = 1.4` and `prob.test = 1..5` via per-test input stubs:
  - `scripts/run_toro1d_uniform_amr.sh` lines 40-43.

## Riemann-problem-based flux construction

- The Euler update uses explicit HLL Riemann flux functions in x/y:
  - `Source/Adv.cpp` lines 98-123 (`hll_flux_x`) and 126-150 (`hll_flux_y`).
- Face fluxes are assembled from left/right reconstructed primitive states and fed to HLL:
  - `Source/Adv.cpp` lines 325-344 (x-faces) and 346-373 (y-faces).

## Second-order reconstruction and time integration

- In non-first-order mode, limited slopes are reconstructed with minmod and positivity-aware scaling:
  - `Source/Adv.cpp` lines 277-323.
- The update is two-stage (predictor-corrector / RK2-like):
  - Stage-1 fluxes from `U^n`: `Source/Adv.cpp` line 377.
  - Stage-2 fluxes from `U*`: `Source/Adv.cpp` line 683.
  - Final average `U^{n+1} = 0.5*(U^n + U^{**})`: `Source/Adv.cpp` lines 704-707 and 721-724.

## Toro Table 4.1 exact-reference generation

- Exact Riemann states for Toro tests 1-5 are hard-coded in postprocessing and used for error computation/overlays:
  - `postprocess_1d.py` lines 249-256.

## Conclusion for marking criteria

- This constitutes a second-order, Riemann-problem-based finite-volume formulation for the 1D Task-1 runs, with exact-solution comparison consistent with the assignment specification.
