# Key Findings (Task 2)

- Task 2 requirement is a 2D verification study: x-split, y-split, and diagonal split on `[0,1]x[0,1]`, each with and without AMR, with clear mesh-patch visualization and one sufficiently accurate resolution.
- We executed the required representative configuration and additionally extended validation to all Toro tests `1..5` in 2D (extra validation; full 2D convergence study not required by the brief).
- Critical correction applied before this final dataset: `ic=1` y-split velocity orientation was fixed by treating Toro states as `(rho, u_n, u_t, p)` and rotating into Cartesian components; this restores proper rotational consistency for nonzero-velocity tests (especially Toro 2 and 5).
- Post-fix rotation-equivalence (`xsplit_vs_ysplit`) is machine-zero for uniform runs and approximately `1e-15` to `1e-13` for AMR runs across tests; see `checks/task2_metrics_by_test.md` and `checks/task2_metrics_long.csv`.
- Diagonal/oblique split results are physically consistent with expected Cartesian grid-alignment diffusion; AMR patch overlays track oblique steep features without visible over-refinement of smooth regions.
- Per-test rho normalization bounds are recorded in `checks/rho_extrema_by_test.tsv`, enabling fair visual comparisons within each test.
- Resolution-matching evidence is explicit: uniform `400x400` vs AMR base `200x200, max_level=1` (effective finest spacing matched); see `checks/resolution_note.txt` and `checks/task2_report_notes.md`.
- This report pack is frozen from the final dataset `results_task2_2d_fullmatrix_20260302_1957` for direct citation in the written assignment.
