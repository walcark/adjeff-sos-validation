"""Everything the comparison asks of adjeff: the depths, and the terms.

Two things come from this side, and both have to be read from adjeff
rather than reproduced:

1. **The optical depths**, out of the Smart-G atmosphere adjeff builds,
   then forced on SOS through ``-AP.MOT`` and ``-AER.AOTref``.  Matching
   the depths by reproducing the Rayleigh formula on both sides would be
   a second place to be wrong; reading one and imposing it on the other
   leaves none.  `depths` also reports ``OD_g``, so the comparison can
   assert there is no gas absorption rather than assume it: Smart-G only
   absorbs through a k-distribution, which adjeff never asks for.
2. **The six 5S terms**, out of `RadiativePipeline`.  One call per sweep
   rather than one per point: the samplers batch over the atmospheric
   dimensions, so a swept AOT costs one Smart-G run per quantity instead
   of one per value.

   Wavelength is the exception.  It is not an atmospheric dimension in
   adjeff but a property of the band, carried by the `ImageDict` key, so
   a wavelength sweep is a sweep over bands: the pipeline is handed all
   of them at once and `run` reassembles the per-band results into a
   ``wl`` dimension.  Passing several wavelengths to `SpectralConfig`
   while handing the pipeline one band computes only that band's,
   silently.

The aerosol optics, the third thing both codes must agree on, are in
:mod:`.optics`.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import xarray as xr
from adjeff.atmosphere import AtmoConfig, GeoConfig, SpectralConfig
from adjeff.atmosphere.atmo_factory import create_atmafgl, grids
from adjeff.core import ImageDict, SensorBand
from adjeff.modules.samplers import RadiativePipeline

from .five_s import QUANTITIES

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
    wl_ref : float
        Reference wavelength of `aot`, in nanometres.  Defaults to
        `WL_REF_NM`, which is adjeff's own.

    Returns
    -------
    Depths
        The depths, read from the atmosphere itself.
    """
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
        afgl_type="afgl_exp_h8km",
        wl_ref=wl_ref,
    )
    ground = atm.sel(z_atm=0.0)

    def scalar_depth(name: str) -> float:
        """Return a ground-level depth as a float, zero when absent."""
        if name not in ground:
            return 0.0
        return float(np.atleast_1d(np.asarray(ground[name].values)).ravel()[0])

    return Depths(
        rayleigh=scalar_depth("OD_r"),
        aerosol=scalar_depth("OD_p"),
        gas=scalar_depth("OD_g"),
        total=scalar_depth("OD_atm"),
        ssa=scalar_depth("ssa_atm"),
    )


def run(
    bands: list[SensorBand],
    aot: float | xr.DataArray,
    sza: float | xr.DataArray,
    vza: float | xr.DataArray,
    raa: float = 180.0,
    rh: float = 50.0,
    href: float = 2.0,
    height: float = 0.0,
    species: dict[str, float] | None = None,
    remove_rayleigh: bool = False,
    n_photons: int | None = None,
) -> xr.Dataset:
    """Run the radiative pipeline and return the six quantities.

    Parameters
    ----------
    bands : list of SensorBand
        Bands to run.  More than one sweeps the wavelength, and the
        result carries a ``wl`` dimension.
    aot : float or xr.DataArray
        Aerosol optical thickness at `WL_REF_NM`.  A ``DataArray`` over
        dim ``"aot"`` sweeps it.
    sza, vza : float or xr.DataArray
        Solar and viewing zenith angles, in degrees.  A ``DataArray``
        over dim ``"sza"`` or ``"vza"`` sweeps them.
    raa : float
        Relative azimuth, in degrees, adjeff's convention: zero is the
        hotspot.
    rh : float
        Relative humidity, in percent.
    href : float
        Aerosol scale height, in km.
    height : float
        Ground elevation, in km.
    species : dict or None
        Aerosol mix.  ``{"sulphate": 1.0}`` when ``None``.
    remove_rayleigh : bool
        Drop molecular scattering.
    n_photons : int or None
        Photons per sampler.  ``None`` keeps adjeff's own defaults, which
        is what makes the comparison a statement about adjeff as shipped.

    Returns
    -------
    xr.Dataset
        The six quantities, over whichever dimensions were swept.
    """
    counts = {}
    if n_photons is not None:
        counts = {
            f"n_ph_{name}": n_photons
            for name in ("sph_alb", "rho_atm", "tdif_up", "tdif_down")
        }

    pipeline = RadiativePipeline(
        atmo_config=AtmoConfig(
            aot=aot,
            rh=rh,
            h=height,
            href=href,
            species=species or {"sulphate": 1.0},
        ),
        geo_config=GeoConfig(sza=sza, vza=vza, saa=raa, vaa=0.0),
        spectral_config=SpectralConfig.from_bands(bands),
        remove_rayleigh=remove_rayleigh,
        **counts,
    )
    scene = pipeline(ImageDict({band: xr.Dataset() for band in bands}))

    per_band = [
        xr.Dataset({name: scene[band][name] for name in QUANTITIES}) for band in bands
    ]
    if len(per_band) == 1:
        return per_band[0]
    wavelengths = [band.wl_nm for band in bands]
    return xr.concat(
        per_band,
        dim=xr.DataArray(wavelengths, dims=["wl"], coords={"wl": wavelengths}),
    )


def scalar(value: xr.DataArray, **where: float) -> float:
    """Select one point of a swept result and return it as a float.

    A coordinate in *where* whose dimension is absent is ignored: a
    quantity that does not depend on a swept parameter does not carry
    its dimension.
    """
    selection = {k: v for (k, v) in where.items() if k in value.dims}
    picked = value.sel(**selection) if selection else value
    flat = np.asarray(picked.values, dtype=float).ravel()
    if flat.size != 1:
        raise ValueError(
            f"selecting {where} on dims {value.dims} left {flat.size} values"
        )
    return float(flat[0])
