r"""The six 5S quantities, from SOS runs at several surface albedos.

SOS prints four of them.  The two it does not, ``rho_atm`` and
``sph_alb``, are read out of the 5S relation itself,

.. math::
    \rho_\mathrm{toa}(A) = \rho_\mathrm{atm}
        + \frac{T^\downarrow T^\uparrow A}{1 - s A}

by running the same atmosphere over several Lambertian albedos ``A``.
A black ground gives ``rho_atm`` directly.  Rearranged for the others,

.. math::
    \frac{\rho_\mathrm{toa}(A) - \rho_\mathrm{atm}}{A}
        = T^\downarrow T^\uparrow + s\,(\rho_\mathrm{toa}(A)
          - \rho_\mathrm{atm})

which is a straight line in ``rho_toa(A) - rho_atm``: its slope is the
spherical albedo and its intercept the two-way transmittance.  Three
albedos make the line overdetermined, and the residual of the fit is
reported.  It is not a free parameter of the comparison but a test of it:
the 5S form is exact over a uniform Lambertian surface, so a residual
above the noise would mean the runs do not describe the atmosphere the
inversion assumes.

The transmittance that comes out is a second, independent reading of
``(tdir_down + tdif_down) * (tdir_up + tdif_up)``, which SOS already
printed term by term.  The two agreeing is what says the inversion is
sound before its ``sph_alb`` is trusted.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .sos import SosCase, SosResult, run

#: The six quantities of the 5S model, under the names adjeff writes.
QUANTITIES = (
    "tdir_down",
    "tdif_down",
    "tdir_up",
    "tdif_up",
    "rho_atm",
    "sph_alb",
)

#: Surface albedos to invert over.  Zero isolates ``rho_atm``; the rest
#: span the range a land scene reaches, and being spread out is what
#: conditions the fit.
ALBEDOS = (0.0, 0.2, 0.5, 0.8)


@dataclass(frozen=True)
class FiveS:
    """The six 5S quantities of one atmosphere and geometry.

    Attributes
    ----------
    tdir_down, tdif_down : float
        Downward transmittances at the solar angle.
    tdir_up, tdif_up : float
        Upward transmittances at the viewing angle.
    rho_atm : float
        Path reflectance over a black ground.
    sph_alb : float
        Spherical albedo of the atmosphere, seen from the ground.
    tau : float
        Total optical depth, as SOS reports it.
    two_way : float
        ``T_down * T_up`` as the inversion reads it.
    two_way_printed : float
        The same product from the four printed transmittances.
    residual : float
        Largest absolute departure of the albedo runs from the fitted
        5S line, in reflectance units.
    sca_ang : float
        Scattering angle SOS states for this geometry, in degrees.
    """

    tdir_down: float
    tdif_down: float
    tdir_up: float
    tdif_up: float
    rho_atm: float
    sph_alb: float
    tau: float
    two_way: float
    two_way_printed: float
    residual: float
    sca_ang: float

    def as_dict(self) -> dict[str, float]:
        """Return the six quantities under the names adjeff writes.

        Returns
        -------
        dict
            Keys ``tdir_down``, ``tdif_down``, ``tdir_up``, ``tdif_up``,
            ``rho_atm`` and ``sph_alb``.
        """
        return {
            "tdir_down": self.tdir_down,
            "tdif_down": self.tdif_down,
            "tdir_up": self.tdir_up,
            "tdif_up": self.tdif_up,
            "rho_atm": self.rho_atm,
            "sph_alb": self.sph_alb,
        }


def five_s(
    case: SosCase,
    vza: float,
    workdir: Path,
    albedos: tuple[float, ...] = ALBEDOS,
) -> FiveS:
    """Run *case* over several albedos and solve for the six quantities.

    Parameters
    ----------
    case : SosCase
        The atmosphere and geometry.  Its ``albedo`` is ignored.
    vza : float
        Viewing zenith angle to report, in degrees.  Must be among
        ``case.vza``.
    workdir : Path
        Directory to run in; one subdirectory per albedo.
    albedos : tuple of float
        Surface albedos, the first of which must be zero.

    Returns
    -------
    FiveS
        The six quantities, plus what it took to get them.

    Raises
    ------
    ValueError
        If *vza* is not among the case's angles, or the albedos do not
        start at zero.
    """
    if vza not in case.vza:
        raise ValueError(f"vza {vza} is not among {case.vza}")
    if albedos[0] != 0.0:
        raise ValueError("the first albedo must be zero, for rho_atm")

    runs: dict[float, SosResult] = {}
    for albedo in albedos:
        tag = f"alb{albedo:.2f}".replace(".", "p")
        case_at = SosCase(**{**case.__dict__, "albedo": albedo})
        runs[albedo] = run(case_at, workdir / tag)

    black = runs[0.0]
    rho_atm = black.rho_toa[vza]

    bright = tuple(a for a in albedos if a != 0.0)
    two_way, sph_alb, residual = invert(
        bright,
        tuple(runs[a].rho_toa[vza] for a in bright),
        rho_atm,
    )

    tdir_up = black.tdir_up(vza)
    tdif_up = black.tdif_up[vza]
    printed = (black.tdir_down + black.tdif_down) * (tdir_up + tdif_up)

    return FiveS(
        tdir_down=black.tdir_down,
        tdif_down=black.tdif_down,
        tdir_up=tdir_up,
        tdif_up=tdif_up,
        rho_atm=rho_atm,
        sph_alb=sph_alb,
        tau=black.tau,
        two_way=two_way,
        two_way_printed=printed,
        residual=residual,
        sca_ang=black.sca_ang[vza],
    )


def invert(
    albedos: tuple[float, ...],
    rho_toa: tuple[float, ...],
    rho_atm: float,
) -> tuple[float, float, float]:
    """Solve the 5S relation for the transmittance and spherical albedo.

    Parameters
    ----------
    albedos : tuple of float
        Surface albedos, all non-zero.
    rho_toa : tuple of float
        Reflectance at the top of the atmosphere, one per albedo.
    rho_atm : float
        Path reflectance, from the black-ground run.

    Returns
    -------
    tuple of float
        ``(two_way, sph_alb, residual)``: the transmittance product, the
        spherical albedo, and the largest departure of the inputs from
        the fitted line, in reflectance units.

    Raises
    ------
    ValueError
        With fewer than two albedos, where the two unknowns are not
        determined.
    """
    if len(albedos) < 2:
        raise ValueError(
            f"two unknowns need at least two bright albedos, got {albedos}"
        )
    bright = np.asarray(albedos, dtype=float)
    excess = np.asarray(rho_toa, dtype=float) - rho_atm
    ratio = excess / bright

    # ratio = two_way + sph_alb * excess, a straight line.
    design = np.column_stack([np.ones_like(excess), excess])
    (two_way, sph_alb), *_ = np.linalg.lstsq(design, ratio, rcond=None)
    residual = float(
        np.max(np.abs(design @ np.array([two_way, sph_alb]) - ratio) * bright)
    )
    return float(two_way), float(sph_alb), residual
