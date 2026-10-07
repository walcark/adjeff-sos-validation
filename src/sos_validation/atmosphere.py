"""The one atmosphere both codes are made to transport light through.

Three things have to agree before a difference between SOS and adjeff
means anything about their methods:

1. **The aerosol optics.**  Handled in :mod:`.optics`: the same CAMS table,
   interpolated the way Smart-G interpolates it.
2. **The optical depths.**  Read here out of the Smart-G atmosphere adjeff
   builds, then forced on SOS through ``-AP.MOT`` and ``-AER.AOTref``.
   Matching the depths by reproducing the Rayleigh formula on both sides
   would be a second place to be wrong; reading one and imposing it on
   the other leaves none.
3. **No gas absorption.**  SOS runs with ``-AP.AbsProfile.Type 7``.
   This one needs nothing on the adjeff side: Smart-G only computes
   gaseous absorption through a k-distribution, which adjeff never asks
   for, so the atmosphere it builds reports ``OD_g`` and ``OD_abs_atm``
   at exactly zero.  `depths` returns both, and the comparison asserts
   they are zero rather than assuming it.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from adjeff.atmosphere.atmo_factory import grids

#: The wavelength adjeff refers its AOT to, from
#: `create_atmosphere(wl_ref_nm=...)`.  Not 550: an AOT of 0.2 means 0.2
#: at 560 nm, so the depth SOS has to be given at another wavelength is
#: not 0.2 either.
WL_REF_NM = 560.0


@dataclass(frozen=True)
class Depths:
    """The optical depths of a Smart-G atmosphere, at the ground.

    Attributes
    ----------
    rayleigh, aerosol, gas, total : float
        ``OD_r``, ``OD_p``, ``OD_g`` and ``OD_atm`` at ``z_atm = 0``.
    ssa : float
        Single scattering albedo of the column.
    """

    rayleigh: float
    aerosol: float
    gas: float
    total: float
    ssa: float


def depths(
    wl_nm: float,
    aot: float,
    rh: float,
    species: dict[str, float],
    href: float = 2.0,
    height: float = 0.0,
    remove_rayleigh: bool = False,
    afgl_type: str = "afgl_exp_h8km",
    wl_ref: float = WL_REF_NM,
) -> Depths:
    """Return the optical depths of the atmosphere adjeff would build.

    The aerosol load adjeff is given is a depth at `wl_ref`, so the depth
    at `wl_nm` is not `aot`: Smart-G rescales it by the extinction ratio.
    SOS has to be told the depth at the simulation wavelength, which is
    what this reads back.

    Parameters
    ----------
    wl_nm : float
        Wavelength, in nanometres.
    aot : float
        Aerosol optical thickness at `wl_ref`.
    rh : float
        Relative humidity, in percent.
    species : dict
        Aerosol mix, proportions summing to one.
    href : float
        Aerosol scale height, in km.
    height : float
        Ground elevation, in km.
    remove_rayleigh : bool
        Drop molecular scattering.
    afgl_type : str
        AFGL profile name.
    wl_ref : float
        Reference wavelength of `aot`, in nanometres.  Defaults to
        `WL_REF_NM`, which is adjeff's own.

    Returns
    -------
    Depths
        The depths, read from the atmosphere itself.
    """
    from adjeff.atmosphere.atmo_factory import create_atmafgl

    # The same vertical grids adjeff's samplers use: they set where the
    # optical depth is sampled, and so its value.
    grid, pfgrid = grids()
    atm = create_atmafgl(
        height=height,
        aot=aot,
        rh=rh,
        wl=wl_nm,
        species=species,
        grid=grid,
        pfgrid=pfgrid,
        zmix=href,
        remove_rayleigh=remove_rayleigh,
        afgl_type=afgl_type,
        wl_ref=wl_ref,
    )
    ground = atm.sel(z_atm=0.0)

    def scalar(name: str) -> float:
        """Return a ground-level depth as a float, zero when absent."""
        if name not in ground:
            return 0.0
        return float(np.atleast_1d(np.asarray(ground[name].values)).ravel()[0])

    return Depths(
        rayleigh=scalar("OD_r"),
        aerosol=scalar("OD_p"),
        gas=scalar("OD_g"),
        total=scalar("OD_atm"),
        ssa=scalar("ssa_atm"),
    )
