# Task 1 Cost Metrics (from existing logs)

| Test | N_eff | Runtime uniform (s) | Runtime AMR (s) | AMR/Uniform runtime | Updated cells uniform | Updated cells AMR | AMR/Uniform updates |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 100 | 0.005966 | 0.006635 | 1.112 | 54800 | 54800 | 1.000 |
| 1 | 200 | 0.018108 | 0.023507 | 1.298 | 220000 | 153072 | 0.696 |
| 1 | 400 | 0.052344 | 0.077873 | 1.488 | 880000 | 717248 | 0.815 |
| 2 | 100 | 0.005483 | 0.005295 | 0.966 | 42000 | 42000 | 1.000 |
| 2 | 200 | 0.013282 | 0.019403 | 1.461 | 166400 | 151952 | 0.913 |
| 2 | 400 | 0.040543 | 0.074138 | 1.829 | 662400 | 914384 | 1.380 |
| 3 | 100 | 0.007301 | 0.007891 | 1.081 | 65200 | 65200 | 1.000 |
| 3 | 200 | 0.020109 | 0.023484 | 1.168 | 258400 | 155072 | 0.600 |
| 3 | 400 | 0.061927 | 0.069852 | 1.128 | 1030400 | 693344 | 0.673 |
| 4 | 100 | 0.007184 | 0.007070 | 0.984 | 60000 | 60000 | 1.000 |
| 4 | 200 | 0.019214 | 0.022140 | 1.152 | 238400 | 146048 | 0.613 |
| 4 | 400 | 0.057757 | 0.063856 | 1.106 | 950400 | 642768 | 0.676 |
| 5 | 100 | 0.011485 | 0.011826 | 1.030 | 106400 | 106400 | 1.000 |
| 5 | 200 | 0.032242 | 0.041061 | 1.274 | 423200 | 351392 | 0.830 |
| 5 | 400 | 0.101910 | 0.148275 | 1.455 | 1684800 | 2010912 | 1.194 |

## Resolution-wise mean ratios

- N_eff=100: mean runtime ratio (AMR/uniform) = **1.034**, mean updated-cell ratio (AMR/uniform) = **1.000**
- N_eff=200: mean runtime ratio (AMR/uniform) = **1.271**, mean updated-cell ratio (AMR/uniform) = **0.730**
- N_eff=400: mean runtime ratio (AMR/uniform) = **1.401**, mean updated-cell ratio (AMR/uniform) = **0.948**

- Interpretation: updated-cell count measures algorithmic work; runtime additionally includes AMR bookkeeping overhead.
