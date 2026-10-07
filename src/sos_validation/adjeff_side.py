"""The same six quantities, from adjeff's radiative samplers.

One `RadiativePipeline` call per sweep rather than one per point: the
samplers batch over the atmospheric dimensions, so a swept AOT costs one
Smart-G run per quantity instead of one per value.

Wavelength is the exception.  It is not an atmospheric dimension in
adjeff but a property of the band, carried by the `ImageDict` key, so a
wavelength sweep is a sweep over bands: the pipeline is handed all of
them at once and this reassembles the per-band results into a ``wl``
dimension.  Passing several wavelengths to `SpectralConfig` while
handing the pipeline one band computes only that band's, silently.
"""

from __future__ import annotations

import numpy as np
import xarray as xr
from adjeff.atmosphere import AtmoConfig, GeoConfig, SpectralConfig
from adjeff.core import ImageDict, SensorBand
from adjeff.modules.samplers import RadiativePipeline

#: The six quantities a forward run produces.
QUANTITIES = (
    "tdir_down",
    "tdif_down",
    "tdir_up",
    "tdif_up",
    "rho_atm",
    "sph_alb",
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
        Aerosol optical thickness at 560 nm.  A ``DataArray`` over dim
        ``"aot"`` sweeps it.
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
    spectral = SpectralConfig.from_bands(bands)

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
        spectral_config=spectral,
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

    Parameters
    ----------
    value : xr.DataArray
        A quantity the pipeline wrote.
    **where : float
        Coordinate values to select, ignored when the dimension is
        absent: a quantity that does not depend on a swept parameter
        does not carry its dimension.

    Returns
    -------
    float
        The selected value.
    """
    selection = {k: v for (k, v) in where.items() if k in value.dims}
    picked = value.sel(**selection) if selection else value
    flat = np.asarray(picked.values, dtype=float).ravel()
    if flat.size != 1:
        raise ValueError(
            f"selecting {where} on dims {value.dims} left {flat.size} values"
        )
    return float(flat[0])
