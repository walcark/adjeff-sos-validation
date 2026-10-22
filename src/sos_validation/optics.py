"""The aerosol optics both codes are given, read from one file.

Smart-G reads a CAMS species from a look-up table under ``SMARTG_DIR_AUXDATA``.
SOS reads an external scattering matrix through ``-AER.Model 4 -AER.ExtData``.
Handing SOS a table of its own (WMO or Shettle & Fenn model) would compare
two aerosol models rather than two radiative transfer methods.

"""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import xarray as xr
from numpy.typing import NDArray

#: Scattering matrix terms, in the order Smart-G stores them along ``stk``.
F11, F21, F33, F34, F22, F44 = range(6)

#: SOS caps an external scattering matrix at 200 angles
# (``CTE_MAXNB_ANG_EXT`` in ``inc/SOS.h``)
SOS_MAX_ANGLES = 200


def auxdata_mixtures() -> Path:
    """Return ``$SMARTG_DIR_AUXDATA/aerosols/OPAC/mixtures``.

    Raises
    ------
    RuntimeError
        If ``SMARTG_DIR_AUXDATA`` is unset, which is also what Smart-G
        itself raises on import.
    """
    root = os.environ.get("SMARTG_DIR_AUXDATA")
    if not root:
        raise RuntimeError(
            "SMARTG_DIR_AUXDATA is unset; it is the directory holding the "
            "CAMS tables that both Smart-G and this comparison read."
        )
    return Path(root) / "aerosols" / "OPAC" / "mixtures"


def sos_angles(n: int = SOS_MAX_ANGLES) -> NDArray[np.float64]:
    """Return the scattering angles the SOS external file is written on.

    Uniform sampling would spend its 200 points where the phase function
    is flat.  ``180 * sin(x)^2`` for ``x`` uniform on ``[0, pi/2]``
    concentrates them at both ends, which is where a coarse mode's
    forward peak and a sphere's backscatter glory live.

    Parameters
    ----------
    n : int
        Number of angles, at most `SOS_MAX_ANGLES`.

    Returns
    -------
    NDArray
        Strictly increasing angles in degrees, from 0 to 180.
    """
    if n > SOS_MAX_ANGLES:
        raise ValueError(f"SOS reads at most {SOS_MAX_ANGLES} angles, got {n}")
    return 180.0 * np.sin(np.linspace(0.0, np.pi / 2, n)) ** 2


@dataclass(frozen=True)
class AerosolOptics:
    """The optics of one CAMS species at one wavelength and humidity.

    Attributes
    ----------
    name : str
        CAMS species, as the table file is named.
    wl_nm : float
        Wavelength, in nanometres.
    rh : float
        Relative humidity, in percent.
    ext : float
        Extinction cross section, in the table's own units.  Only its
        ratio to another wavelength matters: Smart-G scales the profile
        by ``tau_ref``, so any constant factor cancels.
    ssa : float
        Single scattering albedo.
    theta : NDArray
        Scattering angles, in degrees.
    f11, f12, f22, f33 : NDArray
        The four terms SOS reads, on `theta`.
    """

    name: str
    wl_nm: float
    rh: float
    ext: float
    ssa: float
    theta: NDArray[np.float64]
    f11: NDArray[np.float64]
    f12: NDArray[np.float64]
    f22: NDArray[np.float64]
    f33: NDArray[np.float64]

    @classmethod
    def read(
        cls,
        name: str,
        wl_nm: float,
        rh: float,
        n_angles: int = SOS_MAX_ANGLES,
    ) -> AerosolOptics:
        """Interpolate a CAMS table at one wavelength and humidity.

        Parameters
        ----------
        name : str
            CAMS species, e.g. ``"sulphate"``.
        wl_nm : float
            Wavelength in nanometres.
        rh : float
            Relative humidity in percent.
        n_angles : int
            Angles to resample the scattering matrix onto.

        Returns
        -------
        AerosolOptics
            The optics of that species, in SOS's term order.
        """
        path = auxdata_mixtures() / f"{name}_sol.nc"
        with xr.open_dataset(path) as table:
            ds = table.load()

        hum = np.asarray(ds["hum"].values, dtype=float)
        wav = np.asarray(ds["wav"].values, dtype=float)

        def at(values: NDArray[np.float64]) -> NDArray[np.float64]:
            """Interpolate a ``(hum, wav, ...)`` block at (rh, wl_nm)."""
            on_hum = _interp_axis(values, hum, rh, axis=0)
            return _interp_axis(on_hum, wav, wl_nm, axis=0)

        ext = float(at(np.asarray(ds["ext"].values, dtype=float)))
        ssa = float(at(np.asarray(ds["ssa"].values, dtype=float)))
        phase = at(np.asarray(ds["phase"].values, dtype=float))

        theta_in = np.asarray(ds["theta"].values, dtype=float)
        theta = sos_angles(n_angles)

        def resample(term: int, logarithmic: bool) -> NDArray[np.float64]:
            """Move one term from the table's angles onto `theta`."""
            y = phase[term]
            if not logarithmic:
                return np.interp(theta, theta_in, y)
            # interpolating the logarithm of F11 is what keeps the peak from being
            # cut by the coarser grid.
            return np.exp(np.interp(theta, theta_in, np.log(y)))

        terms = {
            "f11": resample(F11, logarithmic=True),
            "f12": resample(F21, logarithmic=False),
            "f22": resample(F22, logarithmic=False),
            "f33": resample(F33, logarithmic=False),
        }
        # Resampling 2000 angles onto 200 shifts the normalisation by a
        # few 1e-4. Scaling all four terms by the same factor restores
        # it without touching the ratios SOS actually reads.
        mu = np.cos(np.deg2rad(theta))
        scale = 1.0 / (0.5 * abs(np.trapezoid(terms["f11"], mu)))

        return cls(
            name=name,
            wl_nm=wl_nm,
            rh=rh,
            ext=ext,
            ssa=ssa,
            theta=theta,
            **{k: v * scale for (k, v) in terms.items()},
        )

    def write_sos_extdata(self, path: Path) -> Path:
        """Write the file ``-AER.ExtData`` reads.

        The format is three labelled scalars, a skipped line, then one row
        per angle: ``theta F11 -F12/F11 F22/F11 F33/F11``.  The sign on
        the second column is SOS's own convention, stated in
        ``SOS_AEROSOLS.F``: it reads a ratio ``-F12/F11``.

        Returns *path*, for chaining.
        """
        path.parent.mkdir(parents=True, exist_ok=True)
        rows = "\n".join(
            f"{t:12.6f} {f11:16.8e} {-f12 / f11:16.8e} "
            f"{f22 / f11:16.8e} {f33 / f11:16.8e}"
            for (t, f11, f12, f22, f33) in zip(
                self.theta, self.f11, self.f12, self.f22, self.f33, strict=True
            )
        )
        path.write_text(
            f"Extinction coefficient : {self.ext:.8e}\n"
            f"Scattering coefficient : {self.ext * self.ssa:.8e}\n"
            f"Number of angles : {self.theta.size}\n"
            "  ANGLE          F11            -F12/F11"
            "         F22/F11          F33/F11\n"
            f"{rows}\n"
        )
        return path


def _interp_axis(
    values: NDArray[np.float64],
    grid: NDArray[np.float64],
    at: float,
    axis: int,
) -> NDArray[np.float64]:
    """Linearly interpolate *values* along *axis*, clamped at both ends.

    *grid* is the increasing coordinate of *axis*, and *axis* is
    collapsed out of the result.
    """
    index = float(np.interp(at, grid, np.arange(grid.size)))
    low = int(np.floor(index))
    high = min(low + 1, grid.size - 1)
    weight = index - low
    moved = np.moveaxis(values, axis, 0)
    return (1.0 - weight) * moved[low] + weight * moved[high]
