# Task 2 Compact Metrics by Test

Metrics are computed from centerline profiles with interpolation to a common coordinate grid.
Rotation check compares x-split against y-split (normal velocity used for y-split).

## Toro test 1

| mode | metric | case | variable | L1 | L2 | Linf |
|---|---|---|---:|---:|---:|---:|
| amr | diag_xy | diag_xslice_vs_yslice | rho | 1.076e-03 | 1.618e-03 | 4.786e-03 |
| amr | exact_1d | xsplit_vs_exact | eint | 1.695e-02 | 6.531e-02 | 5.728e-01 |
| amr | exact_1d | xsplit_vs_exact | p | 2.647e-03 | 7.382e-03 | 1.040e-01 |
| amr | exact_1d | xsplit_vs_exact | rho | 3.977e-03 | 1.052e-02 | 8.093e-02 |
| amr | exact_1d | xsplit_vs_exact | u | 4.946e-03 | 2.766e-02 | 5.749e-01 |
| amr | exact_1d | ysplit_vs_exact | eint | 1.695e-02 | 6.531e-02 | 5.728e-01 |
| amr | exact_1d | ysplit_vs_exact | p | 2.647e-03 | 7.382e-03 | 1.040e-01 |
| amr | exact_1d | ysplit_vs_exact | rho | 3.977e-03 | 1.052e-02 | 8.093e-02 |
| amr | exact_1d | ysplit_vs_exact | u | 4.946e-03 | 2.766e-02 | 5.749e-01 |
| amr | rotation_xy | xsplit_vs_ysplit | eint | 1.100e-15 | 1.828e-15 | 1.465e-14 |
| amr | rotation_xy | xsplit_vs_ysplit | p | 2.121e-16 | 3.623e-16 | 3.525e-15 |
| amr | rotation_xy | xsplit_vs_ysplit | rho | 2.439e-16 | 3.785e-16 | 2.331e-15 |
| amr | rotation_xy | xsplit_vs_ysplit | u | 5.134e-16 | 1.205e-15 | 1.465e-14 |
| uniform | diag_xy | diag_xslice_vs_yslice | rho | 8.468e-04 | 1.668e-03 | 8.296e-03 |
| uniform | exact_1d | xsplit_vs_exact | eint | 1.720e-02 | 6.676e-02 | 7.845e-01 |
| uniform | exact_1d | xsplit_vs_exact | p | 4.677e-03 | 1.393e-02 | 1.921e-01 |
| uniform | exact_1d | xsplit_vs_exact | rho | 5.134e-03 | 1.288e-02 | 1.315e-01 |
| uniform | exact_1d | xsplit_vs_exact | u | 9.517e-03 | 5.131e-02 | 8.485e-01 |
| uniform | exact_1d | ysplit_vs_exact | eint | 1.720e-02 | 6.676e-02 | 7.845e-01 |
| uniform | exact_1d | ysplit_vs_exact | p | 4.677e-03 | 1.393e-02 | 1.921e-01 |
| uniform | exact_1d | ysplit_vs_exact | rho | 5.134e-03 | 1.288e-02 | 1.315e-01 |
| uniform | exact_1d | ysplit_vs_exact | u | 9.517e-03 | 5.131e-02 | 8.485e-01 |
| uniform | rotation_xy | xsplit_vs_ysplit | eint | 0.000e+00 | 0.000e+00 | 0.000e+00 |
| uniform | rotation_xy | xsplit_vs_ysplit | p | 0.000e+00 | 0.000e+00 | 0.000e+00 |
| uniform | rotation_xy | xsplit_vs_ysplit | rho | 0.000e+00 | 0.000e+00 | 0.000e+00 |
| uniform | rotation_xy | xsplit_vs_ysplit | u | 0.000e+00 | 0.000e+00 | 0.000e+00 |

## Toro test 2

| mode | metric | case | variable | L1 | L2 | Linf |
|---|---|---|---:|---:|---:|---:|
| amr | diag_xy | diag_xslice_vs_yslice | rho | 2.994e-03 | 4.108e-03 | 1.380e-02 |
| amr | exact_1d | xsplit_vs_exact | eint | 5.084e-02 | 9.361e-02 | 2.808e-01 |
| amr | exact_1d | xsplit_vs_exact | p | 3.842e-03 | 5.715e-03 | 2.092e-02 |
| amr | exact_1d | xsplit_vs_exact | rho | 8.367e-03 | 1.239e-02 | 3.775e-02 |
| amr | exact_1d | xsplit_vs_exact | u | 1.948e-02 | 2.713e-02 | 7.102e-02 |
| amr | exact_1d | ysplit_vs_exact | eint | 5.084e-02 | 9.361e-02 | 2.808e-01 |
| amr | exact_1d | ysplit_vs_exact | p | 3.842e-03 | 5.715e-03 | 2.092e-02 |
| amr | exact_1d | ysplit_vs_exact | rho | 8.367e-03 | 1.239e-02 | 3.775e-02 |
| amr | exact_1d | ysplit_vs_exact | u | 1.948e-02 | 2.713e-02 | 7.102e-02 |
| amr | rotation_xy | xsplit_vs_ysplit | eint | 1.068e-15 | 1.466e-15 | 5.440e-15 |
| amr | rotation_xy | xsplit_vs_ysplit | p | 3.574e-16 | 6.472e-16 | 3.553e-15 |
| amr | rotation_xy | xsplit_vs_ysplit | rho | 1.061e-15 | 1.708e-15 | 7.772e-15 |
| amr | rotation_xy | xsplit_vs_ysplit | u | 1.195e-15 | 1.507e-15 | 3.997e-15 |
| uniform | diag_xy | diag_xslice_vs_yslice | rho | 1.212e-03 | 1.528e-03 | 4.600e-03 |
| uniform | exact_1d | xsplit_vs_exact | eint | 5.884e-02 | 1.016e-01 | 3.207e-01 |
| uniform | exact_1d | xsplit_vs_exact | p | 6.011e-03 | 9.137e-03 | 2.885e-02 |
| uniform | exact_1d | xsplit_vs_exact | rho | 1.305e-02 | 1.973e-02 | 5.549e-02 |
| uniform | exact_1d | xsplit_vs_exact | u | 3.013e-02 | 3.713e-02 | 9.938e-02 |
| uniform | exact_1d | ysplit_vs_exact | eint | 5.884e-02 | 1.016e-01 | 3.207e-01 |
| uniform | exact_1d | ysplit_vs_exact | p | 6.011e-03 | 9.137e-03 | 2.885e-02 |
| uniform | exact_1d | ysplit_vs_exact | rho | 1.305e-02 | 1.973e-02 | 5.549e-02 |
| uniform | exact_1d | ysplit_vs_exact | u | 3.013e-02 | 3.713e-02 | 9.938e-02 |
| uniform | rotation_xy | xsplit_vs_ysplit | eint | 0.000e+00 | 0.000e+00 | 0.000e+00 |
| uniform | rotation_xy | xsplit_vs_ysplit | p | 0.000e+00 | 0.000e+00 | 0.000e+00 |
| uniform | rotation_xy | xsplit_vs_ysplit | rho | 0.000e+00 | 0.000e+00 | 0.000e+00 |
| uniform | rotation_xy | xsplit_vs_ysplit | u | 0.000e+00 | 0.000e+00 | 0.000e+00 |

## Toro test 3

| mode | metric | case | variable | L1 | L2 | Linf |
|---|---|---|---:|---:|---:|---:|
| amr | diag_xy | diag_xslice_vs_yslice | rho | 4.997e-02 | 2.022e-01 | 1.481e+00 |
| amr | exact_1d | xsplit_vs_exact | eint | 2.958e+01 | 1.541e+02 | 1.512e+03 |
| amr | exact_1d | xsplit_vs_exact | p | 5.034e+00 | 1.443e+01 | 2.773e+02 |
| amr | exact_1d | xsplit_vs_exact | rho | 5.901e-02 | 2.961e-01 | 3.257e+00 |
| amr | exact_1d | xsplit_vs_exact | u | 2.183e-01 | 8.698e-01 | 1.592e+01 |
| amr | exact_1d | ysplit_vs_exact | eint | 2.958e+01 | 1.541e+02 | 1.512e+03 |
| amr | exact_1d | ysplit_vs_exact | p | 5.034e+00 | 1.443e+01 | 2.773e+02 |
| amr | exact_1d | ysplit_vs_exact | rho | 5.901e-02 | 2.961e-01 | 3.257e+00 |
| amr | exact_1d | ysplit_vs_exact | u | 2.183e-01 | 8.698e-01 | 1.592e+01 |
| amr | rotation_xy | xsplit_vs_ysplit | eint | 1.234e-12 | 1.797e-12 | 7.049e-12 |
| amr | rotation_xy | xsplit_vs_ysplit | p | 3.267e-13 | 5.027e-13 | 3.496e-12 |
| amr | rotation_xy | xsplit_vs_ysplit | rho | 1.079e-15 | 3.872e-15 | 3.997e-14 |
| amr | rotation_xy | xsplit_vs_ysplit | u | 1.767e-14 | 2.717e-14 | 2.025e-13 |
| uniform | diag_xy | diag_xslice_vs_yslice | rho | 1.183e-02 | 6.097e-02 | 5.763e-01 |
| uniform | exact_1d | xsplit_vs_exact | eint | 3.474e+01 | 1.775e+02 | 1.608e+03 |
| uniform | exact_1d | xsplit_vs_exact | p | 5.988e+00 | 1.717e+01 | 3.099e+02 |
| uniform | exact_1d | xsplit_vs_exact | rho | 6.107e-02 | 2.817e-01 | 2.935e+00 |
| uniform | exact_1d | xsplit_vs_exact | u | 2.137e-01 | 4.559e-01 | 7.563e+00 |
| uniform | exact_1d | ysplit_vs_exact | eint | 3.474e+01 | 1.775e+02 | 1.608e+03 |
| uniform | exact_1d | ysplit_vs_exact | p | 5.988e+00 | 1.717e+01 | 3.099e+02 |
| uniform | exact_1d | ysplit_vs_exact | rho | 6.107e-02 | 2.817e-01 | 2.935e+00 |
| uniform | exact_1d | ysplit_vs_exact | u | 2.137e-01 | 4.559e-01 | 7.563e+00 |
| uniform | rotation_xy | xsplit_vs_ysplit | eint | 0.000e+00 | 0.000e+00 | 0.000e+00 |
| uniform | rotation_xy | xsplit_vs_ysplit | p | 0.000e+00 | 0.000e+00 | 0.000e+00 |
| uniform | rotation_xy | xsplit_vs_ysplit | rho | 0.000e+00 | 0.000e+00 | 0.000e+00 |
| uniform | rotation_xy | xsplit_vs_ysplit | u | 0.000e+00 | 0.000e+00 | 0.000e+00 |

## Toro test 4

| mode | metric | case | variable | L1 | L2 | Linf |
|---|---|---|---:|---:|---:|---:|
| amr | diag_xy | diag_xslice_vs_yslice | rho | 4.798e-02 | 2.377e-01 | 1.814e+00 |
| amr | exact_1d | xsplit_vs_exact | eint | 2.987e+00 | 1.578e+01 | 1.545e+02 |
| amr | exact_1d | xsplit_vs_exact | p | 4.972e-01 | 1.646e+00 | 3.103e+01 |
| amr | exact_1d | xsplit_vs_exact | rho | 5.885e-02 | 2.995e-01 | 3.569e+00 |
| amr | exact_1d | xsplit_vs_exact | u | 6.984e-02 | 3.052e-01 | 5.230e+00 |
| amr | exact_1d | ysplit_vs_exact | eint | 2.987e+00 | 1.578e+01 | 1.545e+02 |
| amr | exact_1d | ysplit_vs_exact | p | 4.972e-01 | 1.646e+00 | 3.103e+01 |
| amr | exact_1d | ysplit_vs_exact | rho | 5.885e-02 | 2.995e-01 | 3.569e+00 |
| amr | exact_1d | ysplit_vs_exact | u | 6.984e-02 | 3.052e-01 | 5.230e+00 |
| amr | rotation_xy | xsplit_vs_ysplit | eint | 1.107e-13 | 1.695e-13 | 5.969e-13 |
| amr | rotation_xy | xsplit_vs_ysplit | p | 3.260e-14 | 6.728e-14 | 8.313e-13 |
| amr | rotation_xy | xsplit_vs_ysplit | rho | 1.597e-15 | 7.271e-15 | 9.948e-14 |
| amr | rotation_xy | xsplit_vs_ysplit | u | 4.188e-15 | 9.084e-15 | 1.177e-13 |
| uniform | diag_xy | diag_xslice_vs_yslice | rho | 9.486e-03 | 4.603e-02 | 3.659e-01 |
| uniform | exact_1d | xsplit_vs_exact | eint | 3.796e+00 | 1.988e+01 | 1.681e+02 |
| uniform | exact_1d | xsplit_vs_exact | p | 4.791e-01 | 1.511e+00 | 2.778e+01 |
| uniform | exact_1d | xsplit_vs_exact | rho | 5.774e-02 | 2.859e-01 | 3.032e+00 |
| uniform | exact_1d | xsplit_vs_exact | u | 5.262e-02 | 1.218e-01 | 2.024e+00 |
| uniform | exact_1d | ysplit_vs_exact | eint | 3.796e+00 | 1.988e+01 | 1.681e+02 |
| uniform | exact_1d | ysplit_vs_exact | p | 4.791e-01 | 1.511e+00 | 2.778e+01 |
| uniform | exact_1d | ysplit_vs_exact | rho | 5.774e-02 | 2.859e-01 | 3.032e+00 |
| uniform | exact_1d | ysplit_vs_exact | u | 5.262e-02 | 1.218e-01 | 2.024e+00 |
| uniform | rotation_xy | xsplit_vs_ysplit | eint | 0.000e+00 | 0.000e+00 | 0.000e+00 |
| uniform | rotation_xy | xsplit_vs_ysplit | p | 0.000e+00 | 0.000e+00 | 0.000e+00 |
| uniform | rotation_xy | xsplit_vs_ysplit | rho | 0.000e+00 | 0.000e+00 | 0.000e+00 |
| uniform | rotation_xy | xsplit_vs_ysplit | u | 0.000e+00 | 0.000e+00 | 0.000e+00 |

## Toro test 5

| mode | metric | case | variable | L1 | L2 | Linf |
|---|---|---|---:|---:|---:|---:|
| amr | diag_xy | diag_xslice_vs_yslice | rho | 1.263e-01 | 4.606e-01 | 3.626e+00 |
| amr | exact_1d | xsplit_vs_exact | eint | 2.229e+00 | 1.082e+01 | 1.063e+02 |
| amr | exact_1d | xsplit_vs_exact | p | 4.015e+00 | 3.099e+01 | 4.776e+02 |
| amr | exact_1d | xsplit_vs_exact | rho | 2.049e-01 | 9.135e-01 | 8.025e+00 |
| amr | exact_1d | xsplit_vs_exact | u | 3.327e-02 | 3.269e-01 | 7.299e+00 |
| amr | exact_1d | ysplit_vs_exact | eint | 2.229e+00 | 1.082e+01 | 1.063e+02 |
| amr | exact_1d | ysplit_vs_exact | p | 4.015e+00 | 3.099e+01 | 4.776e+02 |
| amr | exact_1d | ysplit_vs_exact | rho | 2.049e-01 | 9.135e-01 | 8.025e+00 |
| amr | exact_1d | ysplit_vs_exact | u | 3.327e-02 | 3.269e-01 | 7.299e+00 |
| amr | rotation_xy | xsplit_vs_ysplit | eint | 9.449e-14 | 2.068e-13 | 1.904e-12 |
| amr | rotation_xy | xsplit_vs_ysplit | p | 8.681e-13 | 1.935e-12 | 1.546e-11 |
| amr | rotation_xy | xsplit_vs_ysplit | rho | 1.268e-14 | 3.516e-14 | 3.730e-13 |
| amr | rotation_xy | xsplit_vs_ysplit | u | 4.369e-15 | 1.135e-14 | 1.597e-13 |
| uniform | diag_xy | diag_xslice_vs_yslice | rho | 1.552e-01 | 4.697e-01 | 2.100e+00 |
| uniform | exact_1d | xsplit_vs_exact | eint | 2.908e+00 | 1.394e+01 | 1.193e+02 |
| uniform | exact_1d | xsplit_vs_exact | p | 1.482e+01 | 1.233e+02 | 1.645e+03 |
| uniform | exact_1d | xsplit_vs_exact | rho | 3.660e-01 | 2.018e+00 | 2.503e+01 |
| uniform | exact_1d | xsplit_vs_exact | u | 1.084e-01 | 1.000e+00 | 1.485e+01 |
| uniform | exact_1d | ysplit_vs_exact | eint | 2.908e+00 | 1.394e+01 | 1.193e+02 |
| uniform | exact_1d | ysplit_vs_exact | p | 1.482e+01 | 1.233e+02 | 1.645e+03 |
| uniform | exact_1d | ysplit_vs_exact | rho | 3.660e-01 | 2.018e+00 | 2.503e+01 |
| uniform | exact_1d | ysplit_vs_exact | u | 1.084e-01 | 1.000e+00 | 1.485e+01 |
| uniform | rotation_xy | xsplit_vs_ysplit | eint | 0.000e+00 | 0.000e+00 | 0.000e+00 |
| uniform | rotation_xy | xsplit_vs_ysplit | p | 0.000e+00 | 0.000e+00 | 0.000e+00 |
| uniform | rotation_xy | xsplit_vs_ysplit | rho | 0.000e+00 | 0.000e+00 | 0.000e+00 |
| uniform | rotation_xy | xsplit_vs_ysplit | u | 0.000e+00 | 0.000e+00 | 0.000e+00 |

