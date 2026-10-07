# adjeff-sos-validation

Validates the six radiative terms of the 5S model that
[adjeff](https://github.com/walcark/Adjeff) computes with the Smart-G
Monte Carlo code against [SOS](https://github.com/CNES/RadiativeTransferCode-SOS),
the CNES successive-orders-of-scattering code MAJA's atmospheric lookup
tables are built with: an independent method, and the reference the
operational processor already uses.

| Term | |
|---|---|
| `tdir_down`, `tdif_down` | Direct and diffuse transmittances, sun to surface |
| `tdir_up`, `tdif_up` | Direct and diffuse transmittances, surface to sensor |
| `rho_atm` | Path reflectance |
| `sph_alb` | Spherical albedo |

## Results

Median relative difference over 339 comparisons: **0.12 %**.

| Term | Worst relative difference |
|---|---|
| `tdir_down`, `tdir_up` | 0.01 % |
| `tdif_down` | 0.79 % |
| `tdif_up` | 0.59 % |
| `sph_alb` | 0.90 % |
| `rho_atm` | 2.0 %, systematically about 0.9 % in a mixed atmosphere |

Five sweeps (AOT, solar angle, viewing angle, wavelength, aerosol scale
height) over three atmospheres (molecules only, aerosol only, both). The
full tables are in [`results/report.md`](results/report.md), every number
in [`results/report.json`](results/report.json).

![adjeff against SOS, one panel per term](results/figures/sos_vs_adjeff.png)

Every point lies on `y = x`. The relative differences, over all 339
comparisons above the 0.01 floor:

![Distribution of the relative differences per term](results/figures/relative_distribution.png)

And against each swept parameter, where the ~1 % on `rho_atm` in a
mixture shows up flat across AOT and angles, and vanishes at a scale
height of 8 km (see [Reading the numbers](#reading-the-numbers)):

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

`plots` needs no GPU: it redraws the figures from any `report.json`
(`--report`, `--out`), in the [SciencePlots](https://github.com/garrettj403/SciencePlots)
`science` style.

`build-sos` clones into `~/.cache/adjeff/sos`, or `SOS_ROOT`. The clone
weighs 4.4 GB, 3.5 GB of which is an absorption database this validation
does not read.

The offline tests need neither SOS, adjeff nor a GPU:

```bash
pixi run -e cpu test
```

## Method

**Same aerosol optics.** SOS reads an external scattering matrix
(`-AER.Model 4`), so both codes read the same CAMS table, interpolated
in wavelength and humidity the way Smart-G's `AerOPAC` does.

**Same optical depths.** The molecular and aerosol depths are read off
the Smart-G atmosphere adjeff builds, then forced on SOS through
`-AP.MOT` and `-AER.AOTref`, rather than recomputed on the SOS side.

**No gas absorption** (`-AP.AbsProfile.Type 7`) and **no phase function
truncation** (`-AER.Tronca 0`), which would renormalise the depth.

**The azimuth.** adjeff's `raa = 0` is the hotspot; SOS's `-SOS.View.Phi 0`
puts sun and sensor in opposite half-planes, so SOS gets `180 - raa`.
Every run re-checks it against the scattering angle SOS prints.

**What SOS gives.** It prints four terms outright. `rho_atm` and `sph_alb`
come from running each atmosphere over four Lambertian albedos and
fitting the 5S relation

```
(rho_toa(A) - rho_atm) / A  =  T_down T_up  +  sph_alb (rho_toa(A) - rho_atm)
```

whose residual (7.6e-07) and two-way transmittance (0.009 % from the
printed one) are reported as checks.

## Reading the numbers

- **adjeff is a Monte Carlo.** At the 2e6 photons used, its noise is
  near 0.1 %.
- **SOS prints its diffuse transmittances to four decimals.** At 1610 nm
  over molecules alone, `tdif_up` is 0.0006 and the rounding alone is
  8 %. Relative differences are ranked only above 0.01; `tdif_down` is
  read from the flux file, one decimal further.
- **The 0.9 % on `rho_atm` in a mixture is the vertical profile.** It is
  not SOS's angular resolution (68 Gauss angles instead of 40 change
  nothing), the forward peak, or the aerosol alone (0.17 %). The two
  codes layer molecules and aerosol on different vertical profiles, and
  SOS takes no keyword to change its own. The scale-height sweep shows
  it: -1.06, -0.93, -0.94, -0.83 % at 0.5, 1, 2 and 4 km, then +0.12 %
  at 8 km, where the aerosol sits with the molecules.

## Layout

| Module | |
|---|---|
| `scripts/build_sos.sh` | Clone and compile SOS |
| `sos_validation/optics.py` | The CAMS optics, and the file SOS reads |
| `sos_validation/atmosphere.py` | The optical depths of adjeff's atmosphere |
| `sos_validation/sos.py` | Run SOS and parse it; the azimuth mapping |
| `sos_validation/five_s.py` | The albedo inversion for `rho_atm` and `sph_alb` |
| `sos_validation/adjeff_side.py` | The six terms from adjeff's `RadiativePipeline` |
| `sos_validation/validate.py` | The sweeps and the command line |
| `sos_validation/report.py` | The Markdown report |
| `sos_validation/plots.py` | The figures |
| `results/` | `report.md`, `report.json` and `figures/` of the last run |

## License

Apache-2.0.
