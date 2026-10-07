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
height of 8 km (see [why](#why-rho_atm-differs-in-a-mixture)):

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

Both codes get the same atmosphere, so only the radiative transfer
can differ:

| Held identical | How |
|---|---|
| Aerosol optics | SOS reads the CAMS scattering matrix Smart-G uses (`-AER.Model 4`), interpolated in wavelength and humidity as `AerOPAC` does |
| Optical depths | Molecular and aerosol depths are read off adjeff's Smart-G atmosphere and forced on SOS (`-AP.MOT`, `-AER.AOTref`) |
| Gas absorption | Off in both (`-AP.AbsProfile.Type 7`) |
| Phase function | Not truncated (`-AER.Tronca 0`), since truncation would renormalise the depth |
| Geometry | SOS gets `180 - raa`: adjeff's `raa = 0` is the hotspot, SOS's is 180. Every run checks this against the scattering angle SOS prints |

SOS prints the four transmittances directly. `rho_atm` and `sph_alb`
come from four SOS runs over Lambertian albedos, fitted to the 5S
relation

```
(rho_toa(A) - rho_atm) / A  =  T_down T_up  +  sph_alb (rho_toa(A) - rho_atm)
```

The fit is checked by its residual (7.6e-07) and by its two-way
transmittance, which is 0.009 % from the one SOS prints.

The score is `(adjeff - SOS) / SOS`, counted only where SOS gives more
than 0.01.

## Limits

- **Monte Carlo noise.** At 2e6 photons, adjeff's noise is about 0.1 %.
- **SOS's printed precision.** SOS prints its diffuse transmittances to
  four decimals. Over molecules alone at 1610 nm, `tdif_up` is 0.0006,
  so rounding alone is worth 8 %; this is why the 0.01 floor exists.
  `tdif_down` is read from the flux file, which has one more decimal.
- **Scope.** The comparison covers Lambertian surfaces, no gas
  absorption, and one CAMS species (sulphate, 50 % humidity): SOS reads
  one scattering matrix, so a mix of species is not compared.

## Why `rho_atm` differs in a mixture

Taken separately, each atmosphere agrees on `rho_atm`: the median
difference is 0.13 % over molecules alone and 0.37 % over aerosol alone.
Mixed together, the difference is a systematic -0.9 %. Three candidate
causes are ruled out:

- SOS's angular resolution: 68 Gauss angles instead of 40 change
  nothing.
- The forward peak of the aerosol phase function: it is not truncated.
- The aerosol optics: they agree when the aerosol is alone.

What remains is the vertical profile. The two codes layer molecules and
aerosol differently, and SOS has no keyword to change its own profile.
`rho_atm` is mostly multiple scattering between the two components, so
it is the most sensitive to which one sits above the other. The
aerosol-scale-height sweep shows the effect:

| Aerosol scale height (km) | 0.5 | 1 | 2 | 4 | 8 |
|---|---|---|---|---|---|
| `rho_atm`, adjeff vs SOS (%) | -1.06 | -0.93 | -0.94 | -0.83 | +0.12 |

The gap closes at 8 km, where the aerosol is spread like the molecules.

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
