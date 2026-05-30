# Task 2 Reporting Notes

## Resolution Matching Statement
- Uniform runs: `amr.n_cell = 400 400`, `amr.max_level = 0`.
- AMR runs: `amr.n_cell = 200 200`, `amr.max_level = 1`, `amr.ref_ratio = 2`.
- Effective finest-grid spacing is matched (`~1/400`) between uniform and AMR.
- `amr.blocking_factor = 8` divides both 400 and 200 exactly.

## Figure Organization for Write-Up
- Main text figures: `figures/main_text/test1/`.
- Appendix figures (additional tests): `figures/appendix/test2/` to `figures/appendix/test5/`.
- Within each test, six figures are provided: `ic=0/1/2`, each with uniform and AMR.

## Metric Files
- Compact markdown tables: `checks/task2_metrics_by_test.md`.
- Full CSV metrics: `checks/task2_metrics_long.csv`.
- Per-test rho normalization: `checks/rho_extrema_by_test.tsv`.
