# Task 3 Quadrant Timing Summary

## Fairness Precautions

- Machine: single laptop node
- Platform: `macOS-26.3-arm64-arm-64bit-Mach-O`
- All timings recorded on the same laptop.
- Quiet timing flags used (`adv.v=0`, `amr.v=0`, `amr.plot_int=1000000`).
- MPI runs used `mpiexec -n p --bind-to none --map-by :OVERSUBSCRIBE`, matching the Task 3 timing scripts and writeup.
- Laptop timings; scaling saturates due to hardware limits/thermals; all comparisons done on same machine.

## Uniform Serial (No AMR)

| N | repeats | mean_s | median_s | std_s |
|---:|---:|---:|---:|---:|
| 768 | 3 | 69.608791 | 69.617823 | 0.152028 |
| 1536 | 3 | 572.865253 | 573.083445 | 0.711657 |
| 3072 | 3 | 4620.244827 | 4595.165286 | 55.256365 |

## Uniform MPI Speedup

| N | cores | repeats | mean_s | median_s | std_s | speedup | efficiency |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 768 | 1 | 3 | 69.608791 | 69.617823 | 0.152028 | 1.000 | 1.000 |
| 768 | 2 | 3 | 38.772523 | 38.882200 | 0.392662 | 1.795 | 0.898 |
| 768 | 4 | 3 | 27.357476 | 27.306983 | 0.951665 | 2.544 | 0.636 |
| 768 | 8 | 3 | 28.418874 | 28.034944 | 1.312129 | 2.449 | 0.306 |
| 768 | 16 | 3 | 34.529022 | 33.823767 | 1.417652 | 2.016 | 0.126 |
| 1536 | 1 | 3 | 572.865253 | 573.083445 | 0.711657 | 1.000 | 1.000 |
| 1536 | 2 | 3 | 345.521426 | 348.823021 | 7.936889 | 1.658 | 0.829 |
| 1536 | 4 | 3 | 255.235778 | 255.272212 | 0.454411 | 2.244 | 0.561 |
| 1536 | 8 | 3 | 218.981420 | 210.921072 | 32.065360 | 2.616 | 0.327 |
| 1536 | 16 | 3 | 206.101579 | 205.988935 | 0.545911 | 2.780 | 0.174 |
| 3072 | 1 | 3 | 4620.244827 | 4595.165286 | 55.256365 | 1.000 | 1.000 |
| 3072 | 2 | 3 | 2566.799285 | 2587.474404 | 38.953993 | 1.800 | 0.900 |
| 3072 | 4 | 3 | 1798.529238 | 1794.469098 | 16.461394 | 2.569 | 0.642 |
| 3072 | 8 | 3 | 1463.201270 | 1473.935160 | 20.448647 | 3.158 | 0.395 |
| 3072 | 16 | 3 | 1522.264132 | 1509.996466 | 21.853377 | 3.035 | 0.190 |

## AMR Match (Serial) Coverage

| effective_N | repeats | mean_s | median_s | std_s | representative L0_cov_% | representative L1_cov_% | representative L2_cov_% |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1536 | 1 | 135.776758 | 135.776758 | 0.000000 | 100.000 | 5.510 | 0.000 |
| 3072 | 3 | 315.041634 | 315.029977 | 0.311533 | 100.000 | 8.160 | 3.470 |

## AMR MPI High-Resolution Speedup

| effective_N | cores | repeats | mean_s | median_s | std_s | speedup | efficiency |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 3072 | 1 | 3 | 315.041634 | 315.029977 | 0.311533 | 1.000 | 1.000 |
| 3072 | 2 | 3 | 181.843167 | 182.757216 | 2.582872 | 1.732 | 0.866 |
| 3072 | 4 | 3 | 138.626330 | 139.066263 | 1.025670 | 2.273 | 0.568 |
| 3072 | 8 | 3 | 140.395056 | 140.763220 | 0.896927 | 2.244 | 0.280 |
| 3072 | 16 | 3 | 166.566436 | 166.966894 | 0.733958 | 1.891 | 0.118 |

## Best High-Resolution AMR vs Uniform

- Uniform N=3072, p=1 mean runtime: **4620.244827 s**
- Fastest measured AMR runtime (effective N=3072): **137.454081 s**
- Fastest measured speed-up vs uniform N=3072 serial: **33.613x**
- Fastest repeated AMR configuration: **p=4**, **138.626330 \pm 1.025670 s** (3 runs)
- Fastest repeated-configuration speed-up vs uniform N=3072 serial: **33.329x**
