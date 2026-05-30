# Full-Field Accuracy Confirmation for Task 3 (Quadrant, t=0.3)

Density is compared against a repeated uniform $N=3072$ serial reference over the full two-dimensional field.
Uniform runs are compared cell-by-cell on their native mesh; AMR runs use valid cells only, with finer levels masking coarse cells, and the uniform reference is bilinearly interpolated to each queried cell centre.

| Family | effective N | L1(rho) | order | L2(rho) | order |
|---|---:|---:|---:|---:|---:|
| Uniform | 768 | 4.345562e-03 | -- | 2.799280e-02 | -- |
| Uniform | 1536 | 2.149105e-03 | 1.016 | 1.614336e-02 | 0.794 |
| AMR | 1536 | 6.129746e-03 | -- | 5.088555e-02 | -- |
| AMR | 3072 | 5.667837e-03 | 0.113 | 5.039762e-02 | 0.014 |

| Matched effective N=1536 comparison | AMR/Uniform L1(rho) ratio |
|---|---:|
| rho | 2.852 |

Interpretation:
- Uniform refinement from 768 to 1536 reduces the full-field density error against the uniform-3072 reference.
- AMR refinement from effective 1536 to effective 3072 also reduces the full-field density error against the same reference.
- At matched effective N=1536, AMR should still be discussed comparatively rather than as uniformly more accurate.
