# Task 1 Key Findings

- Toro Table 4.1 one-dimensional Riemann tests 1-5 were completed with `gamma=1.4`, using three effective resolutions (`N=100,200,400`) for both uniform and AMR runs.
- Uniform-grid validation shows the expected convergence trend: L1 errors in `rho, u, p, e` generally decrease as resolution increases, with test-dependent observed order due to discontinuous solutions.
- Exact solutions were overlaid on all 1D validation plots, and agreement improves visibly with refinement in the dominant wave structures.
- AMR runs used base `N=100` and refinement to matched effective resolutions, then repeated the same error analysis and direct profile comparisons.
- AMR patch locations are now shown directly on AMR-vs-uniform figures using translucent level bands (L1/L2), demonstrating that refinement is concentrated around sharp features and not spread across smooth regions.
- Accuracy at matched high resolution is comparable: across all tests/variables at `N_eff=400`, AMR-to-uniform L1 ratio has mean `0.753` (median `0.729`), with one difficult case (test 4) slightly favoring uniform in some variables.
- Efficiency metrics from logs show nuanced behavior: mean updated-cell ratio (AMR/uniform) is `0.730` at `N_eff=200` and `0.948` at `N_eff=400`, indicating reduced algorithmic work in many cases; however, mean runtime ratio is >1 in this small setup due to AMR overhead.
- Overall, the Task-1 study now supports a high-standard discussion: accuracy validation, convergence quantification, explicit mesh-localization evidence, and transparent cost analysis (work reduction vs wall-clock overhead).
