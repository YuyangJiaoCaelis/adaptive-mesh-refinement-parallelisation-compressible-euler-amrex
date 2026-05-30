# Task 2 Figure/Table Index

This folder is organized for the Task 2 write-up (2D verification). The assignment requires 2D x-split, y-split, and diagonal initial data, each with and without AMR, at one sufficiently accurate resolution. Main-text citations use Toro test 1 as representative; Toro tests 2-5 are provided as appendix validation.

## Main Text (Representative: Toro Test 1)

| Task 2 requirement | Figure(s) to cite |
|---|---|
| 2D x-split (`x<0.5` left, `x>0.5` right), no AMR | `figures/main_text/test1/rho_t1_xsplit_uniform.pdf` (or `.png`) |
| 2D x-split with AMR patch visibility | `figures/main_text/test1/rho_t1_xsplit_amr_with_patches.pdf` (or `.png`) |
| 2D y-split (`y<0.5` left, `y>0.5` right), no AMR | `figures/main_text/test1/rho_t1_ysplit_uniform.pdf` (or `.png`) |
| 2D y-split with AMR patch visibility | `figures/main_text/test1/rho_t1_ysplit_amr_with_patches.pdf` (or `.png`) |
| Oblique/diagonal split, no AMR | `figures/main_text/test1/rho_t1_diag_uniform.pdf` (or `.png`) |
| Oblique/diagonal split with AMR patch visibility | `figures/main_text/test1/rho_t1_diag_amr_with_patches.pdf` (or `.png`) |

## Supporting Evidence (Main Text Tables)

- Resolution matching note: `checks/resolution_note.txt`
- Task-2 report notes: `checks/task2_report_notes.md`
- Compact metrics by test: `checks/task2_metrics_by_test.md`
- Full metrics (CSV): `checks/task2_metrics_long.csv`
- Per-test rho normalization records: `checks/rho_extrema_by_test.tsv`

## Appendix (Additional Validation: Toro Tests 2-5)

- All appendix figures are included under `figures/appendix/test2` ... `figures/appendix/test5`.
- These runs extend Task 2 beyond minimum requirements and provide stronger 2D validation across all Toro tests.
