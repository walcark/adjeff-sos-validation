# Validation of the `adjeff` library with `SOS_ABS_V5.1`

> [!note]
> See [report.md](results/report.md) for the full comparisons.

Validates the six radiative terms of the 5S model that [adjeff](https://github.com/walcark/Adjeff) computes with the Smart-G Monte Carlo code against [SOS](https://github.com/CNES/RadiativeTransferCode-SOS), the CNES successive-orders-of-scattering code MAJA's atmospheric lookup tables are built with: an independent method, and the reference the operational processor already uses.

| Term | |
|---|---|
| `tdir_down`, `tdif_down` | Direct and diffuse transmittances, sun to surface |
| `tdir_up`, `tdif_up` | Direct and diffuse transmittances, surface to sensor |
| `rho_atm` | Path reflectance |
| `sph_alb` | Spherical albedo |

## Results

Median relative difference over 339 comparisons: **0.06 %**.

| Term | Worst relative difference | Where |
|---|---|---|
| `tdir_down`, `tdir_up` | 0.01 % | |
| `rho_atm` | 0.74 % | aerosol only, `sza = 70` |
| `sph_alb` | 0.84 % | aerosol only, 665 nm |
| `tdif_down` | 1.13 % | aerosol only, 1610 nm |
| `tdif_up` | 0.52 % | mixed, 665 nm |

Five sweeps (AOT, solar angle, viewing angle, wavelength, aerosol scale height) over three atmospheres (molecules only, aerosol only, both). The full tables are in [`results/report.md`](results/report.md), every number in [`results/report.json`](results/report.json).

Signed median difference per term and per atmosphere, over the comparisons where SOS gives more than 0.01, in percent:

| Atmosphere | `tdir_down` | `tdif_down` | `tdir_up` | `tdif_up` | `rho_atm` | `sph_alb` | All terms |
|---|---|---|---|---|---|---|---|
| Molecules only | +0.00 | -0.00 | +0.00 | -0.02 | +0.12 | -0.14 | 0.05 |
| Aerosol only | -0.00 | -0.06 | +0.00 | +0.02 | -0.16 | -0.09 | 0.06 |
| Both | -0.00 | -0.02 | +0.00 | +0.10 | -0.21 | -0.08 | 0.06 |

The last column is the median of the absolute difference (in percent) over all six terms. No term is systematically off in any atmosphere: every median is inside the 0.1 % Monte-Carlo noise except `rho_atm`, which is at twice it.

This holds only against an SOS compiled for 500 layers. Against SOS as it ships, the mixed-atmosphere column reads -0.94 % on `rho_atm` and 0.25 % overall, and [that gap is SOS's, not adjeff's](#sos-needs-more-layers-than-it-ships-with).

![adjeff against SOS, one panel per term](results/figures/sos_vs_adjeff.png)

Every point lies on `y = x`. The relative differences, over all 339 comparisons above the 0.01 floor:

![Distribution of the relative differences per term](results/figures/relative_distribution.png)

And against each swept parameter:

![Relative difference against each swept parameter](results/figures/relative_by_sweep.png)

## Running it

Needs a CUDA GPU for adjeff and the Smart-G auxiliary data.

```bash
pixi install
pixi run build-sos                 # clone and compile SOS, once
export SMARTG_DIR_AUXDATA=/path/to/smartg/auxdata
pixi run validate                  # about 30 minutes; writes results/report.*
pixi run validate --sweeps aot --atmospheres aerosol
pixi run plots                     # results/figures/, from results/report.json
```

`plots` needs no GPU: it redraws the figures from any `report.json` (`--report`, `--out`), in the [SciencePlots](https://github.com/garrettj403/SciencePlots) `science` style.

`build-sos` clones into `~/.cache/adjeff/sos`, or `SOS_ROOT`. The clone weighs 4.4 GB, 3.5 GB of which is an absorption database this validation does not read.

Before compiling, it raises `CTE_OS_NT_MIN` in `inc/SOS.h` from the 100 SOS ships with to 500, which is what makes SOS's own vertical discretisation converge (see [below](#sos-needs-more-layers-than-it-ships-with)). `SOS_LAYERS` overrides it; `CTE_OS_NT`, the array bound, allows 600.

The offline tests need neither SOS, adjeff nor a GPU:

```bash
pixi run -e cpu test
```

## Method

Both codes get the same atmosphere, so only the radiative transfer can differ:

| Held identical | How |
|---|---|
| Aerosol optics | SOS reads the CAMS scattering matrix Smart-G uses (`-AER.Model 4`), interpolated in wavelength and humidity as `AerOPAC` does |
| Optical depths | Molecular and aerosol depths are read off adjeff's Smart-G atmosphere and forced on SOS (`-AP.MOT`, `-AER.AOTref`) |
| Gas absorption | Off in both (`-AP.AbsProfile.Type 7`) |
| Phase function | Not truncated (`-AER.Tronca 0`), since truncation would renormalise the depth |
| Geometry | SOS gets `180 - raa`: adjeff's `raa = 0` is the hotspot, SOS's is 180. Every run checks this against the scattering angle SOS prints |
| Vertical resolution | SOS is compiled for 500 layers, not the 100 it ships with, because 101 layers leave its own path reflectance 0.8 % off |

SOS prints the four transmittances directly. `rho_atm` and `sph_alb` come from four SOS runs over Lambertian albedos, fitted to the 5S relation:

```
(rho_toa(A) - rho_atm) / A  =  T_down T_up  +  sph_alb (rho_toa(A) - rho_atm)
```

The fit is checked by its residual (7.6e-07) and by its two-way transmittance, which is 0.009 % from the one SOS prints.

The score is `(adjeff - SOS) / SOS`, counted only where SOS gives more than 0.01.

## Limits

- **Monte Carlo noise.** At 2e6 photons, adjeff's noise is about 0.1 %.
- **SOS's printed precision.** SOS prints its diffuse transmittances to four decimals. Over molecules alone at 1610 nm, `tdif_up` is 0.0006, so rounding alone is worth 8 %; this is why the 0.01 floor exists. `tdif_down` is read from the flux file, which has one more decimal.
- **SOS's vertical resolution is a build-time choice.** SOS as it ships is 0.8 % from its own converged path reflectance in a mixed atmosphere, so `build-sos` recompiles it with 500 layers. Comparing against a stock SOS reproduces that 0.8 % and reads it as an adjeff error; [the section below](#sos-needs-more-layers-than-it-ships-with) is how it was found.
- **Scope.** The comparison covers Lambertian surfaces, no gas absorption, and one CAMS species (sulphate, 50 % humidity): SOS reads one scattering matrix, so a mix of species is not compared.

## SOS needs more layers than it ships with

An earlier run of this comparison, against SOS compiled as it ships, found a systematic -0.9 % on `rho_atm` in a mixed atmosphere and nowhere else: flat across AOT, solar angle, viewing angle and wavelength, absent over molecules alone and over aerosol alone, and vanishing when the aerosol scale height was raised to 8 km. It turned out to be SOS's own vertical discretisation, not a disagreement between the two codes.

SOS cuts the column into layers of equal optical thickness and gives each one a single phase matrix, the Rayleigh and aerosol matrices mixed at that layer's local proportions (`PCAER`, `PCMOL` in its profile log, which run from 0.90 at the ground to 0 at the top for `HA = 2 km`). That is piecewise-constant in a quantity that varies continuously with altitude, so it converges only as the layers get thin. How many layers it uses is a compile-time constant: `CTE_OS_NT_MIN` in `inc/SOS.h`, which gave 101 layers for the tests we made.

Running SOS against with a changing compile-time constant allows to observe the impact of the discretisation:

| `rho_atm` from SOS | 101 layers | 501 layers | Difference |
|---|---|---|---|
| `HA = 2 km`, AOT 0.1 | 0.040075 | 0.039860 | **-0.54 %** |
| `HA = 2 km`, AOT 0.2 | 0.045404 | 0.045046 | **-0.79 %** |
| `HA = 2 km`, AOT 0.4 | 0.057258 | 0.056635 | **-1.09 %** |
| `HA = 8 km`, AOT 0.1 | 0.039916 | 0.039916 | +0.000 % |
| `HA = 8 km`, AOT 0.2 | 0.045229 | 0.045230 | +0.001 % |
| `HA = 8 km`, AOT 0.4 | 0.057102 | 0.057103 | +0.002 % |

At `HA = HR = 8 km` the mixing ratio is the same at every altitude, and therefore the SOS per-layer mixture is exact. This is observed as 101 layers are already converged to 1e-5. As soon as the ratio varies, they are not, and the error is of the size and the sign of the gap that was attributed to adjeff:

| `rho_atm`, `HA = 2 km`, AOT 0.2 | Value | adjeff vs SOS |
|---|---|---|
| adjeff | 0.045008 | |
| SOS, 101 layers | 0.045404 | **-0.87 %** |
| SOS, 501 layers | 0.045046 | **-0.08 %** |

Extrapolating the two layer counts as an error in `1/N` leaves SOS 0.02 % from its own limit at 501 layers, so 0.045046 is the converged value and adjeff was the accurate one all along. Everything above is run against the 501-layer build, and the scale-height sweep that used to show the bias now shows nothing:

## Layout

| Module | |
|---|---|
| `scripts/build_sos.sh` | Clone and compile SOS |
| `sos_validation/optics.py` | The CAMS optics, and the file SOS reads |
| `sos_validation/sos.py` | Run SOS and parse it; the azimuth mapping |
| `sos_validation/five_s.py` | The albedo inversion for `rho_atm` and `sph_alb` |
| `sos_validation/adjeff_side.py` | adjeff's optical depths, and its six terms |
| `sos_validation/validate.py` | The sweeps and the command line |
| `sos_validation/report.py` | The Markdown report |
| `sos_validation/plots.py` | The figures |
| `results/` | `report.md`, `report.json` and `figures/` of the last run |

## License

Apache-2.0.
