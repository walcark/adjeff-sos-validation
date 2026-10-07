"""Checks on the comparison's own arithmetic, without SOS, adjeff or a GPU.

The azimuth mapping and the albedo inversion would shift every number of
the report by a plausible amount rather than break it, so they are
tested against what is known exactly.
"""

from __future__ import annotations

import pytest

from sos_validation.five_s import QUANTITIES, invert
from sos_validation.report import Point, report
from sos_validation.sos import scattering_angle, sos_azimuth

RHO_ATM, TWO_WAY, SPH_ALB = 0.045, 0.875, 0.1235
ALBEDOS = (0.2, 0.5, 0.8)


def forward(albedo: float) -> float:
    """Return the 5S TOA reflectance over a Lambertian *albedo*."""
    return RHO_ATM + TWO_WAY * albedo / (1.0 - SPH_ALB * albedo)


class TestAzimuthConvention:
    """The mapping between SOS's relative azimuth and adjeff's."""

    def test_the_mapping_reproduces_what_sos_printed(self):
        """SOS printed 148.29 degrees at sza 30, vza 1.71, phi 0."""
        assert scattering_angle(30.0, 1.71, 180.0 - 0.0) == pytest.approx(
            148.29, abs=0.01
        )

    def test_the_hotspot_is_at_zero_in_adjeffs_convention(self):
        """``raa = 0`` with ``sza == vza`` is exact backscatter."""
        assert scattering_angle(40.0, 40.0, 0.0) == pytest.approx(180.0)

    def test_sos_and_adjeff_run_backwards_from_each_other(self):
        """SOS's hotspot is at 180, adjeff's at 0."""
        assert sos_azimuth(0.0) == 180.0
        assert sos_azimuth(180.0) == 0.0


class TestInversion:
    """The albedo inversion that gives ``rho_atm`` and ``sph_alb``."""

    def test_it_recovers_the_quantities_it_was_built_from(self):
        """Run the 5S relation forwards, invert it, get the inputs back."""
        rho_toa = tuple(forward(a) for a in ALBEDOS)

        two_way, sph_alb, residual = invert(ALBEDOS, rho_toa, RHO_ATM)

        assert two_way == pytest.approx(TWO_WAY, rel=1e-9)
        assert sph_alb == pytest.approx(SPH_ALB, rel=1e-9)
        assert residual < 1e-12

    def test_the_residual_reports_a_relation_that_does_not_hold(self):
        """A perturbation of 1e-3 on one run has to show in the residual."""
        rho_toa = [forward(a) for a in ALBEDOS]
        rho_toa[1] += 1e-3

        _, _, residual = invert(ALBEDOS, tuple(rho_toa), RHO_ATM)

        assert residual > 1e-4

    def test_one_bright_albedo_is_refused(self):
        """Two unknowns are not determined by one equation."""
        with pytest.raises(ValueError, match="two unknowns"):
            invert((0.5,), (0.5,), RHO_ATM)


def point(value: float, scale: float) -> Point:
    """Return a point where adjeff is *scale* times SOS on every quantity."""
    sos = dict.fromkeys(QUANTITIES, 0.5)
    return Point(
        sweep="aot",
        atmosphere="aerosol",
        varied="aot",
        value=value,
        tau_sos=0.2,
        tau_adjeff=0.2,
        sos=sos,
        adjeff={q: v * scale for q, v in sos.items()},
        sca_ang_sos=150.0,
        sca_ang_geometry=150.0,
        inversion_residual=1e-7,
        two_way_gap=1e-4,
    )


def test_the_report_states_the_relative_differences():
    """A point 1 % above SOS reads ``+1.00`` in its row, and in the median."""
    text = report(
        [point(0.1, 1.01), point(0.2, 1.0)],
        photons=int(2e6),
        gauss=40,
        atmospheres=["aerosol"],
        sweeps=["aot"],
    )

    assert "| 0.1 | +1.00 | +1.00 | +1.00 | +1.00 | +1.00 | +1.00 | 0.2000 |" in text
    assert "## aerosol" in text
    assert "0.50%" in text


def test_the_figures_draw_from_points():
    """The three figures render from report points, with one panel per term."""
    from sos_validation import plots

    points = [point(0.1, 1.01), point(0.2, 1.0)]

    assert len(plots.scatter(points).axes) == len(QUANTITIES)
    assert len(plots.by_sweep(points).axes) == len(QUANTITIES)
    assert len(plots.distribution(points).axes) == 1
