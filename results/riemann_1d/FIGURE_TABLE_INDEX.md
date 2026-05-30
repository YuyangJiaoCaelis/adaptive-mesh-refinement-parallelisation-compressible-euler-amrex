# Task 1 Figure/Table Index

## A. 1D uniform validation (Toro tests 1-5)
Use these figures in the results section:
- `fig_toro1d_test1_uniform.pdf` (or `.png`)
- `fig_toro1d_test2_uniform.pdf` (or `.png`)
- `fig_toro1d_test3_uniform.pdf` (or `.png`)
- `fig_toro1d_test4_uniform.pdf` (or `.png`)
- `fig_toro1d_test5_uniform.pdf` (or `.png`)

What each figure provides:
- 2x2 panels for `rho`, `u`, `p`, `e` vs `x`.
- Uniform `N=100,200,400` on common axes.
- Exact solution (black dashed) overlaid.

Supporting tables for claims:
- `toro1d_errors.md` and `toro1d_errors.csv` for L1 errors and observed orders.

Claims supported:
- L1 errors generally decrease with increasing resolution.
- High-resolution curves are closer to the exact reference.

## B. 1D AMR validation (base N=100)
Use these figures in the AMR section:
- `fig_toro1d_test1_amr_vs_uniform.pdf` (or `.png`)
- `fig_toro1d_test2_amr_vs_uniform.pdf` (or `.png`)
- `fig_toro1d_test3_amr_vs_uniform.pdf` (or `.png`)
- `fig_toro1d_test4_amr_vs_uniform.pdf` (or `.png`)
- `fig_toro1d_test5_amr_vs_uniform.pdf` (or `.png`)

What each figure provides:
- AMR (base=100, refined) vs uniform at matched effective resolution (`N=400`) and exact solution.
- In-panel translucent AMR patch bands (`L1`, `L2`) to make refined regions visible without obscuring solution profiles.
- Bottom annotation of level coverage and refined x-ranges.

Supporting tables for claims:
- `amr_coverage_summary.md` (localization of refinement).
- `toro1d_errors.md` / `toro1d_errors.csv` (accuracy comparison).

Claims supported:
- Refinement remains localized near sharp features; smooth regions remain mostly coarse.
- AMR and uniform solutions are in comparable accuracy range at matched effective resolution.

## C. Efficiency/cost validation (Task 1)
Use these tables:
- `task1_cost_metrics.md` and `task1_cost_metrics.csv`.

What they provide:
- Runtime from solver logs (`Run time = ...`) and total updated-cell counts (`Advanced ... cells`) for uniform vs AMR.
- Ratios by test and resolution (`N_eff=100,200,400`).

Claims supported:
- AMR can reduce algorithmic work (updated-cell count) when refinement stays localized.
- Runtime may still be overhead-dominated for very small 1D-strip jobs; discuss this explicitly.

## D. Numerical-method evidence section
Use:
- `NUMERICAL_METHOD_EVIDENCE.md`

This file provides source-level evidence that the Task-1 runs used a second-order, Riemann-problem-based Euler solver configuration, as required by the brief.
