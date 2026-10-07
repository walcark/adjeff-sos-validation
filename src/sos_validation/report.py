"""The comparison, rendered as Markdown.

Classes
-------
    Point
        One geometry and atmosphere, as both codes answered it.

Functions
---------
    report
        Tables of relative differences, worst cases and self-checks.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .five_s import QUANTITIES
from .sos import TDIF_QUANTUM_FLUX, TDIF_QUANTUM_TRANS

#: Below this, a quantity is too small for a relative difference to mean
#: anything: 1e-4 of Monte-Carlo noise on the molecular diffuse
#: transmittance at 1610 nm, 0.0017, reads as 6 percent.
FLOOR = 0.01

PREAMBLE = """# adjeff against SOS

Relative difference `(adjeff - SOS) / SOS`, in percent, on the six
quantities of the 5S model. Both codes transport light through the same
optical depths and, where there is aerosol, the same scattering matrix.
adjeff ran {photons} photons per Monte-Carlo sampler; SOS projected the
phase function onto {gauss} Gauss angles.

See `README.md` for how the two were made comparable.
"""

WORST_NOTE = """## Worst cases

A relative difference on a quantity near zero says nothing: over
molecules alone at 1610 nm, `tdif_up` is 0.0006, and the rounding step
of the `f7.4` SOS prints it with is 8 percent of that. The relative
column only ranks points where SOS gives more than {floor}; the absolute
column ranks all of them.
"""

MEDIAN_NOTE = """## Summary

Median relative difference over the {count} comparisons above the floor:
{median}%.
"""


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


def _table(header: list[str], rows: list[list[str]]) -> list[str]:
    """Return a Markdown table as lines, with a trailing blank line."""
    return [
        "| " + " | ".join(header) + " |",
        "| " + " | ".join(["---"] * len(header)) + " |",
        *("| " + " | ".join(row) + " |" for row in rows),
        "",
    ]


def _meaningful(point: Point, quantity: str) -> bool:
    """Return whether a relative difference on *quantity* means anything."""
    return bool(np.isfinite(point.relative(quantity))) and (
        abs(point.sos[quantity]) > FLOOR
    )


def report(
    points: list[Point],
    photons: int,
    gauss: int,
    atmospheres: list[str],
    sweeps: list[str],
) -> str:
    """Render *points* as Markdown.

    *photons* and *gauss* are the settings the run used; *atmospheres*
    and *sweeps* set the order of the sections.
    """
    lines = [PREAMBLE.format(photons=f"{photons:.0e}", gauss=gauss)]

    for name in atmospheres:
        mine = [p for p in points if p.atmosphere == name]
        if not mine:
            continue
        lines.append(f"## {name}\n")
        for sweep in sweeps:
            here = [p for p in mine if p.sweep == sweep]
            if not here:
                continue
            lines.append(f"### sweep over {here[0].varied}\n")
            lines += _table(
                [here[0].varied, *QUANTITIES, "tau"],
                [
                    [
                        f"{p.value:g}",
                        *(f"{p.relative(q) * 100:+.2f}" for q in QUANTITIES),
                        f"{p.tau_sos:.4f}",
                    ]
                    for p in here
                ],
            )

    lines.append(WORST_NOTE.format(floor=FLOOR))
    worst = []
    for q in QUANTITIES:
        above = [p for p in points if _meaningful(p, q)]
        rel = max(above, key=lambda p: abs(p.relative(q))) if above else None
        absolute = max(points, key=lambda p: abs(p.absolute(q)))
        worst.append(
            [
                q,
                f"{rel.relative(q) * 100:+.2f}%" if rel else "n/a",
                f"{rel.atmosphere}, {rel.varied}={rel.value:g}" if rel else "n/a",
                f"{absolute.absolute(q):+.2e}",
                f"{absolute.atmosphere}, {absolute.varied}={absolute.value:g}",
            ]
        )
    lines += _table(
        ["quantity", "worst relative", "where", "worst absolute", "where"], worst
    )

    def rounding(quantity: str, quantum: float) -> str:
        """Return the largest share of a value SOS's last decimal is."""
        return f"{100.0 * max(quantum / p.sos[quantity] for p in points):.2f}%"

    lines.append("## Checks on the comparison itself\n")
    lines += _table(
        ["check", "worst"],
        [
            [
                "optical depth gap between the two codes",
                f"{max(abs(p.tau_sos - p.tau_adjeff) for p in points):.2e}",
            ],
            [
                "departure of the SOS albedo runs from the 5S line",
                f"{max(p.inversion_residual for p in points):.2e}",
            ],
            [
                "inverted against printed two-way transmittance",
                f"{max(abs(p.two_way_gap) for p in points) * 100:.3f}%",
            ],
            [
                "azimuth convention, SOS scattering angle against the geometry",
                f"{max(abs(p.sca_ang_sos - p.sca_ang_geometry) for p in points):.2f}"
                " deg",
            ],
            [
                "worst `tdif_up` rounding, SOS printing it `f7.4`",
                rounding("tdif_up", TDIF_QUANTUM_TRANS),
            ],
            [
                "worst `tdif_down` rounding, the flux file `f9.5`",
                rounding("tdif_down", TDIF_QUANTUM_FLUX),
            ],
        ],
    )

    above_floor = [
        abs(p.relative(q)) for p in points for q in QUANTITIES if _meaningful(p, q)
    ]
    lines.append(
        MEDIAN_NOTE.format(
            count=len(above_floor), median=f"{np.median(above_floor) * 100:.2f}"
        )
    )
    return "\n".join(lines)
