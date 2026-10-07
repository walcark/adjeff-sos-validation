# adjeff against SOS

Relative difference `(adjeff - SOS) / SOS`, in percent, on the six
quantities of the 5S model. Both codes transport light through the same
optical depths and, where there is aerosol, the same scattering matrix.
adjeff ran 2e+06 photons per Monte-Carlo sampler; SOS projected the
phase function onto 40 Gauss angles.

See `README.md` for how the two were made comparable.

## rayleigh

### sweep over sza

| sza | tdir_down | tdif_down | tdir_up | tdif_up | rho_atm | sph_alb | tau |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | -0.00 | -0.14 | +0.00 | +0.23 | +0.58 | -0.59 | 0.0902 |
| 20 | +0.00 | +0.21 | +0.00 | +0.23 | +0.39 | -0.60 | 0.0902 |
| 40 | -0.00 | +0.06 | +0.00 | +0.23 | -0.49 | -0.60 | 0.0902 |
| 60 | -0.00 | -0.18 | +0.00 | +0.23 | +0.02 | -0.60 | 0.0902 |
| 70 | -0.01 | +0.01 | +0.00 | +0.23 | +0.47 | -0.59 | 0.0902 |

### sweep over vza

| vza | tdir_down | tdif_down | tdir_up | tdif_up | rho_atm | sph_alb | tau |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | +0.00 | -0.05 | +0.00 | +0.12 | -0.13 | -0.37 | 0.0902 |
| 10 | +0.00 | -0.05 | +0.00 | -0.21 | -0.29 | -0.38 | 0.0902 |
| 20 | +0.00 | -0.05 | +0.00 | +0.12 | +0.12 | -0.37 | 0.0902 |
| 40 | +0.00 | -0.05 | +0.00 | +0.05 | +0.09 | -0.37 | 0.0902 |
| 60 | +0.00 | -0.05 | +0.01 | -0.27 | -0.04 | -0.37 | 0.0902 |

### sweep over wl_nm

| wl_nm | tdir_down | tdif_down | tdir_up | tdif_up | rho_atm | sph_alb | tau |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 443 | +0.00 | -0.04 | -0.00 | -0.08 | -0.10 | -0.18 | 0.2358 |
| 560 | +0.00 | -0.49 | +0.00 | +0.24 | +0.02 | +0.08 | 0.0902 |
| 665 | -0.00 | +0.56 | -0.00 | -0.03 | +0.19 | -0.01 | 0.0448 |
| 865 | -0.00 | -1.29 | +0.00 | -0.37 | +0.27 | -0.15 | 0.0155 |
| 1610 | +0.00 | +0.30 | +0.00 | +7.72 | -1.25 | +1.43 | 0.0013 |

## aerosol

### sweep over aot

| aot | tdir_down | tdif_down | tdir_up | tdif_up | rho_atm | sph_alb | tau |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 0.05 | +0.00 | -0.07 | +0.00 | -0.18 | +0.53 | -0.69 | 0.0500 |
| 0.1 | +0.01 | -0.28 | +0.00 | +0.31 | +0.28 | -0.02 | 0.1000 |
| 0.2 | -0.00 | -0.06 | +0.00 | +0.07 | +0.81 | -0.00 | 0.2000 |
| 0.4 | -0.00 | -0.05 | +0.00 | +0.08 | -0.30 | -0.33 | 0.4000 |
| 0.8 | +0.01 | -0.01 | +0.00 | -0.05 | +0.02 | -0.30 | 0.8000 |

### sweep over sza

| sza | tdir_down | tdif_down | tdir_up | tdif_up | rho_atm | sph_alb | tau |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | +0.00 | -0.00 | +0.00 | +0.11 | +0.09 | -0.60 | 0.2000 |
| 20 | -0.00 | -0.28 | +0.00 | +0.11 | +0.49 | -0.60 | 0.2000 |
| 40 | +0.00 | -0.08 | +0.00 | +0.11 | +0.41 | -0.60 | 0.2000 |
| 60 | +0.00 | +0.04 | +0.00 | +0.11 | +0.37 | -0.60 | 0.2000 |
| 70 | +0.01 | +0.09 | +0.00 | +0.11 | -0.63 | -0.59 | 0.2000 |

### sweep over vza

| vza | tdir_down | tdif_down | tdir_up | tdif_up | rho_atm | sph_alb | tau |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | -0.00 | -0.28 | +0.00 | -0.07 | +0.18 | -0.25 | 0.2000 |
| 10 | -0.00 | -0.28 | +0.00 | -0.03 | +0.67 | -0.25 | 0.2000 |
| 20 | -0.00 | -0.28 | +0.00 | -0.26 | +0.29 | -0.25 | 0.2000 |
| 40 | -0.00 | -0.28 | +0.00 | +0.15 | -0.48 | -0.25 | 0.2000 |
| 60 | -0.00 | -0.28 | +0.00 | -0.09 | +0.52 | -0.25 | 0.2000 |

### sweep over wl_nm

| wl_nm | tdir_down | tdif_down | tdir_up | tdif_up | rho_atm | sph_alb | tau |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 443 | -0.00 | +0.18 | +0.00 | +0.40 | -0.13 | -0.19 | 0.2572 |
| 560 | -0.00 | +0.13 | +0.00 | -0.16 | +0.22 | +0.26 | 0.2000 |
| 665 | +0.01 | +0.16 | -0.00 | +0.04 | -0.33 | -0.21 | 0.1558 |
| 865 | -0.00 | -0.40 | -0.00 | +0.12 | -0.51 | -0.79 | 0.0968 |
| 1610 | -0.00 | -0.27 | +0.00 | +0.59 | +0.43 | +0.58 | 0.0203 |

## both

### sweep over aot

| aot | tdir_down | tdif_down | tdir_up | tdif_up | rho_atm | sph_alb | tau |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 0.05 | -0.00 | +0.49 | +0.00 | +0.53 | -0.44 | +0.14 | 0.1402 |
| 0.1 | -0.01 | +0.13 | +0.00 | +0.25 | -0.65 | -0.25 | 0.1902 |
| 0.2 | -0.00 | +0.27 | +0.00 | +0.21 | -1.02 | -0.78 | 0.2902 |
| 0.4 | -0.00 | +0.30 | +0.00 | +0.33 | -1.28 | -0.66 | 0.4902 |
| 0.8 | -0.01 | +0.15 | +0.00 | +0.12 | -1.12 | -0.49 | 0.8902 |

### sweep over sza

| sza | tdir_down | tdif_down | tdir_up | tdif_up | rho_atm | sph_alb | tau |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | +0.01 | +0.16 | +0.00 | +0.47 | -0.89 | -0.64 | 0.2902 |
| 20 | +0.00 | +0.21 | +0.00 | +0.47 | -0.99 | -0.64 | 0.2902 |
| 40 | -0.00 | +0.13 | +0.00 | +0.47 | -1.15 | -0.64 | 0.2902 |
| 60 | +0.00 | +0.21 | +0.00 | +0.47 | -0.28 | -0.64 | 0.2902 |
| 70 | +0.00 | +0.12 | +0.00 | +0.47 | -0.37 | -0.64 | 0.2902 |

### sweep over vza

| vza | tdir_down | tdif_down | tdir_up | tdif_up | rho_atm | sph_alb | tau |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | -0.00 | +0.29 | +0.00 | +0.30 | -0.76 | -0.71 | 0.2902 |
| 10 | -0.00 | +0.29 | +0.00 | +0.21 | -1.10 | -0.70 | 0.2902 |
| 20 | -0.00 | +0.29 | +0.00 | +0.08 | -0.78 | -0.70 | 0.2902 |
| 40 | -0.00 | +0.29 | +0.00 | +0.37 | -0.78 | -0.70 | 0.2902 |
| 60 | -0.00 | +0.29 | +0.01 | +0.02 | -0.43 | -0.71 | 0.2902 |

### sweep over wl_nm

| wl_nm | tdir_down | tdif_down | tdir_up | tdif_up | rho_atm | sph_alb | tau |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 443 | +0.00 | +0.25 | -0.00 | +0.46 | -0.57 | -0.23 | 0.4930 |
| 560 | -0.00 | +0.11 | +0.00 | +0.11 | -1.04 | -0.65 | 0.2902 |
| 665 | +0.00 | +0.21 | -0.00 | -0.17 | -1.33 | -0.26 | 0.2006 |
| 865 | +0.01 | +0.79 | -0.00 | +0.38 | -2.02 | -0.90 | 0.1123 |
| 1610 | +0.00 | +0.04 | -0.00 | +0.39 | +0.06 | -0.28 | 0.0215 |

### sweep over href

| href | tdir_down | tdif_down | tdir_up | tdif_up | rho_atm | sph_alb | tau |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 0.5 | -0.00 | +0.22 | +0.00 | +0.32 | -1.26 | -0.52 | 0.2902 |
| 1 | -0.00 | +0.41 | +0.00 | +0.39 | -1.15 | -0.73 | 0.2902 |
| 2 | -0.00 | +0.30 | +0.00 | +0.27 | -0.87 | -0.59 | 0.2902 |
| 4 | -0.00 | +0.37 | +0.00 | +0.16 | -1.04 | -0.46 | 0.2902 |
| 8 | -0.00 | -0.26 | +0.00 | +0.06 | +0.16 | +0.03 | 0.2902 |

## Worst cases

A relative difference on a quantity near zero says nothing: over
molecules alone at 1610 nm, `tdif_up` is 0.0006, and the rounding step
of the `f7.4` SOS prints it with is 8 percent of that. The relative
column only ranks points where SOS gives more than 0.01; the absolute
column ranks all of them.

| quantity | worst relative | where | worst absolute | where |
| --- | --- | --- | --- | --- |
| tdir_down | -0.01% | both, aot=0.8 | +4.97e-05 | both, wl_nm=865 |
| tdif_down | +0.79% | both, wl_nm=865 | +1.03e-03 | both, aot=0.4 |
| tdir_up | +0.01% | both, vza=60 | +5.64e-05 | rayleigh, vza=60 |
| tdif_up | +0.59% | aerosol, wl_nm=1610 | +1.21e-03 | both, wl_nm=443 |
| rho_atm | -2.02% | both, wl_nm=865 | -9.52e-04 | both, aot=0.8 |
| sph_alb | -0.90% | both, wl_nm=865 | -1.09e-03 | both, aot=0.8 |

## Checks on the comparison itself

| check | worst |
| --- | --- |
| optical depth gap between the two codes | 4.72e-05 |
| departure of the SOS albedo runs from the 5S line | 7.55e-07 |
| inverted against printed two-way transmittance | 0.009% |
| azimuth convention, SOS scattering angle against the geometry | 0.00 deg |
| worst `tdif_up` rounding, SOS printing it `f7.4` | 8.33% |
| worst `tdif_down` rounding, the flux file `f9.5` | 0.68% |

## Summary

Median relative difference over the 339 comparisons above the floor:
0.12%.
