# Accuracy Confirmation for Task 3 (Quadrant, t=0.3)

Reference solution for this centreline comparison: **Repeated uniform reference (N=3072, p=1)**.
A directly post-processed uniform $N=3072$ serial plotfile is used as the high-resolution centreline reference for this comparison.

| Family | effective N | L1 error: rho(x, y=0.5) | L2 error: rho(x, y=0.5) | L1 error: rho(x=0.5, y) | L2 error: rho(x=0.5, y) |
|---|---:|---:|---:|---:|---:|
| AMR | 1536 | 6.041530e-03 | 1.848722e-02 | 6.506777e-03 | 2.255616e-02 |
| AMR | 3072 | 4.217700e-03 | 2.721798e-02 | 4.708724e-03 | 2.564903e-02 |
| Uniform | 768 | 4.340423e-03 | 1.732041e-02 | 4.340423e-03 | 1.732041e-02 |
| Uniform | 1536 | 3.329921e-03 | 1.283560e-02 | 3.329921e-03 | 1.283560e-02 |

| Family | Resolution step | L1 ratio x | L1 order x | L1 ratio y | L1 order y |
|---|---|---:|---:|---:|---:|
| Uniform | 768 $\rightarrow$ 1536 | 1.303 | 0.382 | 1.303 | 0.382 |
| AMR | 1536 $\rightarrow$ 3072 | 1.432 | 0.518 | 1.382 | 0.467 |

| Matched effective N=1536 comparison | AMR/Uniform L1 ratio x | AMR/Uniform L1 ratio y |
|---|---:|---:|
| rho | 1.814 | 1.954 |

Interpretation:
- Uniform refinement from 768 to 1536 reduces the centerline density L1 error in both diagnostics.
- AMR refinement from effective 1536 to effective 3072 also reduces the direct centerline error against the high-resolution uniform reference.
- At matched effective N=1536, AMR is not uniformly better on every centerline metric, so the Task 3 accuracy statement should remain comparative rather than absolute.
