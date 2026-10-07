"""Compare adjeff's radiative samplers against SOS, parameter by parameter.

Five sweeps, each moving one thing away from a reference point:

=========  ====================================================
``aot``    the aerosol load, at nadir
``sza``    the solar angle, which ``tdir_down`` and
           ``tdif_down`` follow
``vza``    the viewing angle, which ``tdir_up`` and ``tdif_up``
           follow
``wl``     the wavelength, over the Sentinel-2 range
``href``   the aerosol scale height, which only ``rho_atm``
           in a mixture is sensitive to
=========  ====================================================

Each is run over three atmospheres, because a disagreement in the mixed
one says nothing about where it comes from:

===============  ===================================================
``rayleigh``     molecules only, no phase function to interpolate
``aerosol``      aerosol only, testing the scattering matrix both
                 codes were handed
``both``         what a user actually runs
===============  ===================================================

Usage
-----
    pixi run validate
    pixi run validate --sweeps aot --atmospheres aerosol
"""

from __future__ import annotations

import argparse
import json
import tempfile
import time
from dataclasses import asdict, dataclass, field, replace
from pathlib import Path

import xarray as xr
from adjeff.core import S2Band

from . import adjeff_side, atmosphere, sos
from .adjeff_side import QUANTITIES
from .five_s import FiveS, five_s
from .optics import AerosolOptics
from .report import report

#: Sentinel-2 bands the wavelength sweep runs over.  In adjeff a
#: wavelength is a property of the band, so a wavelength sweep is a band
#: sweep.
BANDS = (S2Band.B01, S2Band.B03, S2Band.B04, S2Band.B8A, S2Band.B11)

#: Band of every other sweep, 560 nm.
BASE_BAND = S2Band.B03


@dataclass(frozen=True)
class Params:
    """One point of the parameter space, angles in degrees."""

    wl_nm: float = BASE_BAND.wl_nm
    aot: float = 0.2
    sza: float = 30.0
    vza: float = 0.0
    raa: float = 180.0
    rh: float = 50.0
    href: float = 2.0
    species: dict[str, float] = field(default_factory=lambda: {"sulphate": 1.0})


#: Reference point every sweep passes through.
BASE = Params()

#: Sweep name to ``(varied parameter, values)``.
SWEEPS: dict[str, tuple[str, tuple[float, ...]]] = {
    "aot": ("aot", (0.05, 0.1, 0.2, 0.4, 0.8)),
    "sza": ("sza", (0.0, 20.0, 40.0, 60.0, 70.0)),
    "vza": ("vza", (0.0, 10.0, 20.0, 40.0, 60.0)),
    "wl": ("wl_nm", tuple(band.wl_nm for band in BANDS)),
    "href": ("href", (0.5, 1.0, 2.0, 4.0, 8.0)),
}

#: Dimension each varied parameter lands on in the adjeff result.
DIMS = {"wl_nm": "wl", "aot": "aot", "sza": "sza", "vza": "vza", "href": "href"}

#: Atmosphere name to ``(molecules, aerosol)``.
ATMOSPHERES: dict[str, tuple[bool, bool]] = {
    "rayleigh": (True, False),
    "aerosol": (False, True),
    "both": (True, True),
}

#: Sweeps an atmosphere does not depend on.
SKIP: set[tuple[str, str]] = {
    ("rayleigh", "aot"),
    ("rayleigh", "href"),
    ("aerosol", "href"),
}

#: Photons per Monte-Carlo sampler: a tenth of adjeff's default, which
#: puts the noise near 0.1 percent and a full run under an hour.
N_PHOTONS = int(2e6)


@dataclass
class Point:
    """One geometry and atmosphere, as both codes answered it.

    Attributes
    ----------
    sweep, atmosphere : str
        Keys of `SWEEPS` and `ATMOSPHERES`.
    varied : str
        Name of the varied parameter.
    value : float
        Its value at this point.
    tau_sos, tau_adjeff : float
        Total optical depth each code used; forced equal.
    sos, adjeff : dict
        The six quantities.
    sca_ang_sos, sca_ang_geometry : float
        Scattering angle SOS states, and the one the geometry implies in
        adjeff's convention: their agreement checks the azimuth mapping.
    inversion_residual : float
        Departure of the SOS albedo runs from the 5S line.
    two_way_gap : float
        Relative gap between the inverted and the printed two-way
        transmittance.
    """

    sweep: str
    atmosphere: str
    varied: str
    value: float
    tau_sos: float
    tau_adjeff: float
    sos: dict[str, float]
    adjeff: dict[str, float]
    sca_ang_sos: float
    sca_ang_geometry: float
    inversion_residual: float
    two_way_gap: float

    def absolute(self, quantity: str) -> float:
        """Return ``adjeff - sos`` on *quantity*."""
        return self.adjeff[quantity] - self.sos[quantity]

    def relative(self, quantity: str) -> float:
        """Return ``(adjeff - sos) / sos`` on *quantity*, ``nan`` if SOS is 0."""
        reference = self.sos[quantity]
        if reference == 0.0:
            return float("nan")
        return (self.adjeff[quantity] - reference) / reference


def sos_point(
    params: Params,
    molecules: bool,
    aerosol: bool,
    depths: atmosphere.Depths,
    workdir: Path,
) -> FiveS:
    """Run SOS at *params*, through the optical depths adjeff computed.

    Raises
    ------
    ValueError
        For an aerosol mix of more than one species: SOS reads one
        external scattering matrix.
    """
    optics = None
    if aerosol:
        if len(params.species) != 1:
            raise ValueError(
                "SOS reads one external scattering matrix, so the comparison "
                f"takes one species at a time, got {sorted(params.species)}"
            )
        optics = AerosolOptics.read(next(iter(params.species)), params.wl_nm, params.rh)

    case = sos.SosCase(
        wl_nm=params.wl_nm,
        aot=depths.aerosol if aerosol else 0.0,
        mot=depths.rayleigh if molecules else 0.0,
        sza=params.sza,
        aerosol=optics,
        vza=(params.vza,),
        raa=sos.sos_azimuth(params.raa),
        aer_scale_height=params.href,
    )
    return five_s(case, params.vza, workdir)


def adjeff_points(sweep: str, molecules: bool, aerosol: bool) -> xr.Dataset:
    """Run adjeff once over a whole sweep.

    ``aerosol=False`` is expressed as a negligible AOT, which is what
    adjeff offers.

    Raises
    ------
    ValueError
        For an atmosphere with neither molecules nor aerosol.
    """
    if not molecules and not aerosol:
        raise ValueError("an atmosphere with neither molecules nor aerosol")

    varied, values = SWEEPS[sweep]
    swept: dict[str, xr.DataArray] = {}
    if sweep != "wl":
        dim = DIMS[varied]
        swept[varied] = xr.DataArray(
            list(values), dims=[dim], coords={dim: list(values)}
        )

    def pick(name: str, default: float) -> float | xr.DataArray:
        return swept.get(name, default)

    return adjeff_side.run(
        bands=list(BANDS) if sweep == "wl" else [BASE_BAND],
        aot=pick("aot", BASE.aot if aerosol else 1e-9),
        sza=pick("sza", BASE.sza),
        vza=pick("vza", BASE.vza),
        raa=BASE.raa,
        rh=BASE.rh,
        href=pick("href", BASE.href),  # type: ignore[arg-type]
        species=BASE.species,
        remove_rayleigh=not molecules,
        n_photons=N_PHOTONS,
    )


def collect(sweeps: list[str], atmospheres: list[str], workdir: Path) -> list[Point]:
    """Run both codes over the requested sweeps and atmospheres.

    Returns
    -------
    list of Point
        One per (atmosphere, sweep, value), in that order.
    """
    points: list[Point] = []
    for name in atmospheres:
        molecules, aerosol = ATMOSPHERES[name]
        for sweep in sweeps:
            if (name, sweep) in SKIP:
                print(f"skip   {name:9s} {sweep:4s} (not a parameter of it)")
                continue
            varied, values = SWEEPS[sweep]
            started = time.monotonic()
            batched = adjeff_points(sweep, molecules, aerosol)
            print(
                f"adjeff {name:9s} {sweep:4s} {time.monotonic() - started:6.1f}s",
                flush=True,
            )

            for value in values:
                params = replace(BASE, **{varied: value})
                depths = atmosphere.depths(
                    wl_nm=params.wl_nm,
                    aot=params.aot if aerosol else 0.0,
                    rh=params.rh,
                    species=params.species,
                    href=params.href,
                    remove_rayleigh=not molecules,
                )
                if depths.gas != 0.0:
                    raise RuntimeError(
                        f"the Smart-G atmosphere absorbs (OD_g={depths.gas}), "
                        "which SOS was told not to"
                    )

                tag = f"{name}-{sweep}-{value:g}".replace(".", "p")
                ref = sos_point(params, molecules, aerosol, depths, workdir / tag)
                got = {
                    q: adjeff_side.scalar(batched[q], **{DIMS[varied]: value})
                    for q in QUANTITIES
                }
                point = Point(
                    sweep=sweep,
                    atmosphere=name,
                    varied=varied,
                    value=float(value),
                    tau_sos=ref.tau,
                    tau_adjeff=depths.total,
                    sos=ref.as_dict(),
                    adjeff=got,
                    sca_ang_sos=ref.sca_ang,
                    sca_ang_geometry=sos.scattering_angle(
                        params.sza, params.vza, params.raa
                    ),
                    inversion_residual=ref.residual,
                    two_way_gap=ref.two_way / ref.two_way_printed - 1.0,
                )
                points.append(point)
                print(
                    f"  {sweep}={value:<7g} "
                    + "  ".join(
                        f"{q}={point.relative(q) * 100:+6.2f}%" for q in QUANTITIES
                    ),
                    flush=True,
                )
    return points


def main() -> None:
    """Parse the command line, run the comparison, write the report."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--sweeps", nargs="+", choices=sorted(SWEEPS), default=list(SWEEPS)
    )
    parser.add_argument(
        "--atmospheres",
        nargs="+",
        choices=sorted(ATMOSPHERES),
        default=list(ATMOSPHERES),
    )
    parser.add_argument(
        "--gauss",
        type=int,
        default=sos.N_GAUSS_RAD,
        help=(
            f"Gauss angles SOS projects the phase function onto, at most "
            f"{sos.MAX_GAUSS_RAD}; two runs at two values separate SOS's "
            "angular resolution from a real disagreement."
        ),
    )
    parser.add_argument(
        "--workdir",
        type=Path,
        default=Path(tempfile.gettempdir()) / "adjeff-sos",
        help="where the SOS runs happen; wiped per case",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("results"),
        help="directory report.md and report.json are written to",
    )
    args = parser.parse_args()
    sos.set_gauss(args.gauss, min(args.gauss, sos.MAX_GAUSS_AER))

    points = collect(args.sweeps, args.atmospheres, args.workdir)
    args.out.mkdir(parents=True, exist_ok=True)
    raw, md = args.out / "report.json", args.out / "report.md"
    raw.write_text(json.dumps([asdict(p) for p in points], indent=2) + "\n")
    md.write_text(
        report(
            points,
            photons=N_PHOTONS,
            gauss=sos.N_GAUSS_RAD,
            atmospheres=list(ATMOSPHERES),
            sweeps=list(SWEEPS),
        )
    )
    print(f"\nwrote {md} and {raw}")


if __name__ == "__main__":
    main()
