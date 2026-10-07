"""Drive the SOS executable and read back what it computed.

SOS is run once per (atmosphere, geometry, surface albedo).  Four of the
six 5S quantities come straight out of its ``-SOS.Trans`` file:

===================  =========================================
5S quantity          Where SOS states it
===================  =========================================
``tdir_down``        ``Direct transmission TOA -> surface``
``tdif_down``        ``Downward diffuse flux at BOA``, in the flux
                     file.  The transmission file states the same
                     number as ``td(thetas)`` but to four decimals
                     against five, and `TDIF_QUANTUM` is what that
                     costs.
``tdif_up``          ``td(thetav)``, one row per angle, four decimals
``tdir_up``          analytic, from the profile's total depth
===================  =========================================

The other two, ``rho_atm`` and ``sph_alb``, are not printed and come from
the albedo inversion in :mod:`.five_s`.

Everything runs without gas absorption (``-AP.AbsProfile.Type 7``) and
without phase function truncation (``-AER.Tronca 0``).  Truncation would
move the optical depth SOS transports light through, so ``tdir`` would no
longer be ``exp(-tau/mu)`` for the ``tau`` the comparison declares.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from .optics import AerosolOptics

#: Where `scripts/build_sos.sh` puts the clone and the executable.
SOS_ROOT = Path(os.environ.get("SOS_ROOT", Path.home() / ".cache/adjeff/sos"))

#: Gauss angles for the radiance field and for the phase matrix.  The
#: phase matrix is capped at 100 by ``CTE_MIE_NBMU_MAX``.  The radiance
#: field looks capped at 80 by ``CTE_OS_NBMU_MAX``, but the real limit is
#: 68: ``SOS_ANGLES`` derives ``OS_NS = 2 * NbGauss`` and checks it
#: against ``CTE_OS_NS_MAX``, which ships at 136.  This is the
#: resolution SOS projects the phase function onto, so it sets the floor
#: on how well a forward peak can agree with Smart-G's 2000-point table.
#: `set_gauss` moves it, which is how the comparison tells its own
#: resolution from a real disagreement.
N_GAUSS_RAD = 40
N_GAUSS_AER = 40

#: What the shipped ``inc/SOS.h`` actually allows.
MAX_GAUSS_RAD = 68
MAX_GAUSS_AER = 100


def set_gauss(radiance: int, aerosol: int) -> None:
    """Change the Gauss angle counts for subsequent runs.

    Parameters
    ----------
    radiance : int
        Angles for the radiance field, at most `MAX_GAUSS_RAD`.
    aerosol : int
        Angles for the phase matrix, at most 100.

    Raises
    ------
    ValueError
        Beyond what SOS was compiled for, where it would overrun its
        own tables rather than refuse.
    """
    global N_GAUSS_RAD, N_GAUSS_AER
    if not 1 <= radiance <= MAX_GAUSS_RAD or not 1 <= aerosol <= MAX_GAUSS_AER:
        raise ValueError(
            f"SOS allows at most {MAX_GAUSS_RAD} radiance and "
            f"{MAX_GAUSS_AER} aerosol Gauss angles, got {radiance} and "
            f"{aerosol}; beyond that it stops on ERROR_4000 rather than "
            "overrunning its tables"
        )
    N_GAUSS_RAD, N_GAUSS_AER = radiance, aerosol


#: Rounding step of the diffuse transmittances SOS prints.  The
#: transmission file writes them ``f7.4`` and the flux file ``f9.5``, so
#: a value of 0.0007 carries one significant digit and a relative
#: difference on it is meaningless.  Reported alongside, so a percentage
#: can be read against what it is measured on.
TDIF_QUANTUM_TRANS = 0.5e-4
TDIF_QUANTUM_FLUX = 0.5e-5

#: Maximal scattering order.  The default is 100; at an optical depth
#: below one the series converges in far fewer, and the cost of asking
#: for more is nil.
N_ORDERS = 100


class SosError(RuntimeError):
    """SOS ran but did not produce the files it was asked for."""


@dataclass(frozen=True)
class SosCase:
    """One SOS simulation.

    Attributes
    ----------
    wl_nm : float
        Wavelength, in nanometres.
    aot : float
        Aerosol optical thickness at `wl_nm`.  SOS is given the reference
        wavelength equal to the simulation wavelength, which
        ``-AER.Model 4`` requires.
    mot : float
        Molecular optical thickness, forced through ``-AP.MOT`` so that
        both codes transport light through the same Rayleigh depth.  Zero
        isolates the aerosol.
    sza : float
        Solar zenith angle, in degrees.
    aerosol : AerosolOptics or None
        The species.  ``None`` runs a purely molecular atmosphere.
    albedo : float
        Lambertian surface albedo.
    vza : tuple of float
        Viewing zenith angles to report, added to the Gauss set with a
        zero weight so they change no integral.
    raa : float
        Relative azimuth, in SOS's convention: 180 degrees puts the
        satellite and the sun in the same half-plane.
    aer_scale_height : float
        Aerosol scale height, in km.
    mol_scale_height : float
        Molecular scale height, in km.  8 km is what
        ``afgl_exp_h8km`` uses on the adjeff side.
    truncate : bool
        Let SOS truncate the forward peak of the phase function.  Off by
        default: truncation renormalises the optical depth with the peak,
        so ``tdir`` would no longer be ``exp(-tau/mu)`` for the ``tau``
        the comparison declares.  Turning it on is how the Legendre
        reconstruction of an untruncated peak is told from a real
        disagreement on ``rho_atm``.
    """

    wl_nm: float
    aot: float
    mot: float
    sza: float
    aerosol: AerosolOptics | None = None
    albedo: float = 0.0
    vza: tuple[float, ...] = (0.0,)
    raa: float = 180.0
    aer_scale_height: float = 2.0
    mol_scale_height: float = 8.0
    truncate: bool = False


@dataclass
class SosResult:
    """What one SOS run says about the atmosphere it simulated.

    Attributes
    ----------
    tau : float
        Total optical depth of the profile, molecular plus aerosol.
    tdir_down : float
        Direct transmittance from TOA to the surface, at the solar angle.
    tdif_down : float
        Diffuse downward transmittance, at the solar angle.
    tdif_up : dict
        Diffuse upward transmittance, keyed by viewing zenith angle.
    rho_toa : dict
        Reflectance at the top of the atmosphere, keyed by viewing
        zenith angle, in Smart-G's convention ``pi L / (Esun mu_s)``.
    sca_ang : dict
        Scattering angle SOS states for each viewing angle, in degrees.
        Kept so the azimuth convention is checked rather than assumed.
    flux_down_direct, flux_down_diffuse : float
        Downward fluxes at the surface, normalised by the TOA solar flux.
    flux_up_diffuse : float
        Upward diffuse flux at TOA, same normalisation.
    """

    tau: float
    tdir_down: float
    tdif_down: float
    tdif_up: dict[float, float] = field(default_factory=dict)
    rho_toa: dict[float, float] = field(default_factory=dict)
    sca_ang: dict[float, float] = field(default_factory=dict)
    flux_down_direct: float = float("nan")
    flux_down_diffuse: float = float("nan")
    flux_up_diffuse: float = float("nan")

    def tdir_up(self, vza: float) -> float:
        """Return the direct upward transmittance at *vza*.

        SOS does not print it, and it needs no model: a beam leaving the
        surface at *vza* is attenuated by ``exp(-tau / cos(vza))``, with
        the same ``tau`` SOS reports for the profile.

        Parameters
        ----------
        vza : float
            Viewing zenith angle, in degrees.

        Returns
        -------
        float
            The direct upward transmittance.
        """
        return float(np.exp(-self.tau / np.cos(np.deg2rad(vza))))


def executable() -> Path:
    """Return the SOS executable, raising if it has not been built.

    Returns
    -------
    Path
        Path to ``SOS_ABS_MAIN.exe``.
    """
    exe = SOS_ROOT / "exe" / "SOS_ABS_MAIN.exe"
    if not exe.is_file():
        raise SosError(
            f"{exe} is missing; run `pixi run build-sos` (or set SOS_ROOT)."
        )
    return exe


def run(case: SosCase, workdir: Path) -> SosResult:
    """Run one SOS simulation and parse its outputs.

    Parameters
    ----------
    case : SosCase
        The simulation to run.
    workdir : Path
        Directory to run in.  Wiped first: SOS appends to some of its
        files and reuses a binary intermediate, so a stale directory
        silently mixes two runs.

    Returns
    -------
    SosResult
        The parsed quantities.

    Raises
    ------
    SosError
        If SOS fails, or writes no transmission file.
    """
    if workdir.exists():
        shutil.rmtree(workdir)
    (workdir / "MIE").mkdir(parents=True)
    (workdir / "SURF").mkdir(parents=True)

    angles = workdir / "user_angles.txt"
    angles.write_text("".join(f"{v:.6f}\n" for v in sorted(case.vza)))

    argv = [
        str(executable()),
        "-SOS_Main.Wa",
        f"{case.wl_nm / 1000.0:.6f}",
        "-SOS_Main.ResRoot",
        str(workdir),
        "-SOS_Main.Log",
        "main.log",
        "-ANG.Rad.NbGauss",
        str(N_GAUSS_RAD),
        "-ANG.Aer.NbGauss",
        str(N_GAUSS_AER),
        "-ANG.Rad.UserAngFile",
        str(angles),
        "-ANG.Thetas",
        f"{case.sza:.6f}",
        "-ANG.Log",
        "angles.log",
        "-SOS.View",
        "1",
        "-SOS.View.Phi",
        f"{case.raa:.6f}",
        "-SOS.IGmax",
        str(N_ORDERS),
        "-SOS.Ipolar",
        "1",
        "-SOS.ResFileUp",
        "up.txt",
        "-SOS.ResFileDown",
        "down.txt",
        "-SOS.Trans",
        "trans.txt",
        "-SOS.Flux",
        "flux.txt",
        "-SOS.Log",
        "sos.log",
        # No gas absorption: the two codes are being compared on how they
        # transport scattered light, and their gas profiles are unrelated.
        "-AP.AbsProfile.Type",
        "7",
        "-AP.MOT",
        f"{case.mot:.8f}",
        "-AP.HR",
        f"{case.mol_scale_height:.4f}",
        "-AP.AerProfile.Type",
        "1",
        "-AP.AerHS.HA",
        f"{case.aer_scale_height:.4f}",
        "-AP.Log",
        "profile.log",
        "-SURF.Dir",
        str(workdir / "SURF"),
        "-SURF.Type",
        "0",
        "-SURF.Alb",
        f"{case.albedo:.6f}",
        "-SURF.Log",
        "surf.log",
    ]

    if case.aerosol is None:
        # SOS-ABS carries no -AP.AOT: the load is set by -AER.AOTref
        # alone, and zero there is a purely molecular atmosphere.
        argv += ["-AER.AOTref", "0.0"]
    else:
        extdata = case.aerosol.write_sos_extdata(workdir / "extdata.txt")
        argv += [
            "-AER.DirMie",
            str(workdir / "MIE"),
            "-AER.Model",
            "4",
            "-AER.ExtData",
            str(extdata),
            "-AER.Waref",
            f"{case.wl_nm / 1000.0:.6f}",
            "-AER.AOTref",
            f"{case.aot:.8f}",
            # Truncation renormalises the phase function and the optical
            # depth with it, which would break tdir = exp(-tau/mu).
            "-AER.Tronca",
            "0",
            "-AER.ResFile",
            "aer.txt",
            "-AER.Log",
            "aer.log",
        ]

    # SOS reads SOS_ABS_ROOT itself, to locate the data it ships with.
    # Without it SOS_AEROSOLS stops at ERROR_925, having already written
    # a plausible-looking log.
    env = dict(os.environ, SOS_ABS_ROOT=str(SOS_ROOT))
    proc = subprocess.run(argv, capture_output=True, text=True, cwd=workdir, env=env)
    trans = workdir / "SOS" / "trans.txt"
    if not trans.is_file():
        raise SosError(
            f"SOS wrote no transmission file.\n"
            f"stdout tail:\n{proc.stdout[-2000:]}\n"
            f"stderr tail:\n{proc.stderr[-2000:]}"
        )

    result = _parse_trans(trans.read_text(), case.vza)
    _parse_flux((workdir / "SOS" / "flux.txt").read_text(), result)
    # The flux file states the same downward diffuse transmittance to one
    # more decimal, so it replaces what the transmission file printed.
    result.tdif_down = result.flux_down_diffuse
    _parse_radiance((workdir / "SOS" / "up.txt").read_text(), case, result)
    return result


def _parse_trans(text: str, vza: tuple[float, ...]) -> SosResult:
    """Read the direct and diffuse transmittances SOS prints.

    Parameters
    ----------
    text : str
        Contents of the ``-SOS.Trans`` file.
    vza : tuple of float
        Angles to pick out of the upward table.

    Returns
    -------
    SosResult
        Filled except for the flux and radiance fields.
    """
    direct = _one_float(text, r"Direct transmission TOA -> surface\s*:\s*([-\d.E+]+)")
    td_s = _one_float(text, r"thetas\s*=\s*[-\d.]+\s+td\(thetas\)\s*=\s*([-\d.E+]+)")

    table = {
        float(a): float(t)
        for (a, t) in re.findall(
            r"thetav\s*=\s*([-\d.]+)\s+td\(thetav\)\s*=\s*([-\d.E+]+)", text
        )
    }
    result = SosResult(tau=float("nan"), tdir_down=direct, tdif_down=td_s)
    for angle in vza:
        result.tdif_up[angle] = _nearest(table, angle, "td(thetav)")
    return result


def _parse_flux(text: str, result: SosResult) -> None:
    """Read the fluxes and the total optical depth into *result*.

    Parameters
    ----------
    text : str
        Contents of the ``-SOS.Flux`` file.
    result : SosResult
        Updated in place.
    """
    result.flux_down_direct = _one_float(
        text, r"Downward direct flux at BOA\s*:\s*([-\d.E+]+)"
    )
    result.flux_down_diffuse = _one_float(
        text, r"Downward diffuse flux at BOA\s*:\s*([-\d.E+]+)"
    )
    result.flux_up_diffuse = _one_float(
        text, r"Upward diffuse flux at TOA[^:]*:\s*([-\d.E+]+)"
    )
    # The profile table repeats one row per layer; every row carries the
    # column-integrated depth, so the first is enough.
    rows = re.findall(
        r"^\s*[-\d.]+\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s*$",
        text,
        re.MULTILINE,
    )
    if not rows:
        raise SosError("no profile table in the flux file")
    result.tau = float(rows[0][3])


def _parse_radiance(text: str, case: SosCase, result: SosResult) -> None:
    """Read the TOA reflectance at the requested viewing angles.

    SOS normalises its Stokes I by the solar irradiance alone,
    ``pi L / Esun``, where Smart-G divides by ``Esun * cos(sza)`` as
    well.  The conversion is that cosine.

    The sign of the printed VZA selects the half-plane: a positive angle
    is the azimuth given by ``-SOS.View.Phi``, a negative one that
    azimuth plus 180 degrees.  This reads the positive side.

    Parameters
    ----------
    text : str
        Contents of the ``-SOS.ResFileUp`` file.
    case : SosCase
        The case, for the solar angle.
    result : SosResult
        Updated in place.
    """
    rows: dict[float, float] = {}
    scattering: dict[float, float] = {}
    for line in text.splitlines():
        if line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 3:
            continue
        angle = float(parts[0])
        if angle < 0.0:
            continue
        rows[angle] = float(parts[2])
        scattering[angle] = float(parts[1])

    mu_s = np.cos(np.deg2rad(case.sza))
    for angle in case.vza:
        result.rho_toa[angle] = _nearest(rows, angle, "radiance") / mu_s
        result.sca_ang[angle] = _nearest(scattering, angle, "scattering angle")


def _one_float(text: str, pattern: str) -> float:
    """Return the single float *pattern* captures in *text*.

    Parameters
    ----------
    text : str
        Text to search.
    pattern : str
        Regular expression with one capturing group.

    Returns
    -------
    float
        The captured value.

    Raises
    ------
    SosError
        If the pattern does not match.
    """
    found = re.search(pattern, text)
    if found is None:
        raise SosError(f"{pattern!r} not found in the SOS output")
    return float(found.group(1))


def _nearest(table: dict[float, float], angle: float, what: str) -> float:
    """Return *table* at the key closest to *angle*, within 0.01 degrees.

    The tolerance is what turns a silently wrong angle into a failure:
    a user angle SOS did not honour would otherwise be answered with its
    nearest Gauss node, degrees away.

    Parameters
    ----------
    table : dict
        Angle to value.
    angle : float
        Angle wanted, in degrees.
    what : str
        Name of the quantity, for the error message.

    Returns
    -------
    float
        The value at *angle*.

    Raises
    ------
    SosError
        If no key is within 0.01 degrees.
    """
    if not table:
        raise SosError(f"no {what} rows in the SOS output")
    key = min(table, key=lambda k: abs(k - angle))
    if abs(key - angle) > 0.01:
        raise SosError(
            f"{what} is tabulated at {key} degrees, not the {angle} asked "
            "for; check that the user angle file reached SOS"
        )
    return table[key]


def scattering_angle(sza: float, vza: float, raa: float) -> float:
    """Return the scattering angle of a geometry, in adjeff's convention.

    ``raa = 0`` is the hotspot: at ``sza == vza`` it gives 180 degrees.
    adjeff's test suite settled its sensor azimuth against this relation,
    in ``test_the_path_reflectance_follows_the_scattering_angle``.

    Parameters
    ----------
    sza, vza, raa : float
        Solar zenith, viewing zenith and relative azimuth, in degrees.

    Returns
    -------
    float
        The scattering angle, in degrees.
    """
    s, v, r = np.deg2rad([sza, vza, raa])
    cos_theta = -np.cos(s) * np.cos(v) - np.sin(s) * np.sin(v) * np.cos(r)
    return float(np.rad2deg(np.arccos(np.clip(cos_theta, -1.0, 1.0))))


def sos_azimuth(raa: float) -> float:
    """Convert a relative azimuth from adjeff's convention to SOS's.

    SOS states its own: ``0`` degrees puts satellite and sun in opposite
    half-planes, which is adjeff's ``raa = 180``.  The two run backwards
    from each other, so ``-SOS.View.Phi`` takes ``180 - raa``.

    Measured rather than taken on trust: at ``sza = 30``, ``vza = 1.71``
    and ``-SOS.View.Phi 0``, SOS prints a scattering angle of 148.29
    degrees, and `scattering_angle` gives 148.3 for ``raa = 180``.
    `SosResult.sca_ang` carries SOS's value so every run re-checks it.

    Parameters
    ----------
    raa : float
        Relative azimuth in adjeff's convention, in degrees.

    Returns
    -------
    float
        The value to pass as ``-SOS.View.Phi``.
    """
    return 180.0 - raa
