"""Figures of the comparison, drawn from ``report.json``.

Functions
---------
    scatter
        adjeff against SOS, one panel per quantity.
    by_sweep
        Relative difference against each swept parameter.
    distribution
        Relative differences of each quantity, over every point.
    main
        Command line: draw the three figures.

Usage
-----
    pixi run plots
    pixi run plots --report results/report.json --out results/figures
"""

from __future__ import annotations

import argparse
import functools
import json
from collections.abc import Callable
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import scienceplots  # noqa: E402, F401  (registers the styles)
from matplotlib.figure import Figure  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

from .five_s import QUANTITIES  # noqa: E402
from .report import FLOOR, Point  # noqa: E402

#: The SciencePlots sheets.  ``no-latex`` spares a TeX install; drop it
#: for TeX-rendered text.
SHEETS = ["science", "no-latex"]

#: Atmosphere to ``(label, colour, marker)``, colours from SciencePlots'
#: ``vibrant``.  Colour and shape both carry identity, so the figures read
#: in greyscale and for any vision; markers are hollow, so overlapping
#: points of different atmospheres stay visible.
STYLE = {
    "rayleigh": ("Molecules only", "#0077BB", "o"),
    "aerosol": ("Aerosol only", "#EE7733", "s"),
    "both": ("Molecules and aerosol", "#EE3377", "D"),
}

#: Sweep to its axis label.
SWEEP_LABELS = {
    "aot": "AOT at 560 nm",
    "sza": "Solar zenith angle (°)",
    "vza": "Viewing zenith angle (°)",
    "wl": "Wavelength (nm)",
    "href": "Aerosol scale height (km)",
}

#: The relative-difference axis of `by_sweep`, in percent; beyond it a
#: point is drawn on the edge, as a triangle pointing out.  No
#: atmosphere uses the triangle, so the two cannot be confused.
LIMIT = 3.0

INK = "#0b0b0b"
MUTED = "#52514e"
BAND = "#ebebeb"


def load(path: Path) -> list[Point]:
    """Return the points of a ``report.json``."""
    return [Point(**entry) for entry in json.loads(path.read_text())]


def _styled(draw: Callable[[list[Point]], Figure]) -> Callable[[list[Point]], Figure]:
    """Draw *draw*'s figure under the SciencePlots `SHEETS`."""

    @functools.wraps(draw)
    def wrapper(points: list[Point]) -> Figure:
        with plt.style.context(SHEETS):
            return draw(points)

    return wrapper


def _meaningful(point: Point, quantity: str) -> bool:
    """Return whether SOS is far enough from zero for a relative difference."""
    return abs(point.sos[quantity]) > FLOOR


def _hollow(atmosphere: str, size: float) -> dict:
    """Return the keywords of an atmosphere's hollow marker for `Axes.plot`."""
    _, colour, marker = STYLE[atmosphere]
    return {
        "linestyle": "",
        "marker": marker,
        "markersize": size,
        "markerfacecolor": "none",
        "markeredgecolor": colour,
        "markeredgewidth": 1.3,
    }


def _legend(fig: Figure, atmospheres: list[str], filled: bool = False) -> None:
    """Add the atmosphere legend above the panels, with the figure's markers."""
    handles = []
    for a in atmospheres:
        marker = _hollow(a, 7)
        if filled:
            marker["markerfacecolor"] = marker["markeredgecolor"]
        handles.append(Line2D([], [], label=STYLE[a][0], **marker))
    fig.legend(
        handles=handles,
        loc="lower center",
        ncol=len(handles),
        frameon=False,
        bbox_to_anchor=(0.5, 1.0),
        fontsize=10,
    )


@_styled
def scatter(points: list[Point]) -> Figure:
    """Return adjeff against SOS, one panel per quantity, with ``y = x``.

    Each atmosphere gets its own marker size, largest first, so markers
    that overlap on ``y = x`` nest instead of hiding each other.
    """
    atmospheres = [a for a in STYLE if any(p.atmosphere == a for p in points)]
    sizes = dict(zip(atmospheres, (9.5, 6.5, 4.0), strict=False))
    fig, axes = plt.subplots(2, 3, figsize=(10, 6.8))
    for ax, q in zip(axes.flat, QUANTITIES, strict=True):
        sos = np.array([p.sos[q] for p in points])
        adj = np.array([p.adjeff[q] for p in points])
        low, high = min(sos.min(), adj.min()), max(sos.max(), adj.max())
        pad = 0.05 * (high - low)
        edges = [low - pad, high + pad]
        ax.plot(edges, edges, color=MUTED, lw=0.8, ls="--", zorder=1)
        for a in atmospheres:
            mine = [p for p in points if p.atmosphere == a]
            ax.plot(
                [p.sos[q] for p in mine],
                [p.adjeff[q] for p in mine],
                zorder=3,
                **_hollow(a, sizes[a]),
            )
        ax.set_xlim(edges)
        ax.set_ylim(edges)
        ax.set_aspect("equal")
        ax.set_title(q, fontsize=11)
        ax.set_xlabel("SOS")
        ax.set_ylabel("adjeff")

        rel = [abs(p.relative(q)) for p in points if _meaningful(p, q)]
        if rel:
            ax.text(
                0.96,
                0.05,
                f"median {np.median(rel) * 100:.2f} %\nworst {max(rel) * 100:.2f} %",
                transform=ax.transAxes,
                ha="right",
                va="bottom",
                color=MUTED,
                fontsize=8.5,
            )
    fig.tight_layout()
    _legend(fig, atmospheres)
    return fig


@_styled
def by_sweep(points: list[Point]) -> Figure:
    """Return the relative difference against each swept parameter.

    Rows are quantities, columns sweeps.  Faded markers are points where
    SOS gives less than `FLOOR`, whose relative difference means little.
    """
    sweeps = [s for s in SWEEP_LABELS if any(p.sweep == s for p in points)]
    atmospheres = [a for a in STYLE if any(p.atmosphere == a for p in points)]
    fig, axes = plt.subplots(
        len(QUANTITIES),
        len(sweeps),
        figsize=(2.7 * len(sweeps), 1.7 * len(QUANTITIES)),
        sharey=True,
        sharex="col",
        squeeze=False,
    )
    for row, q in enumerate(QUANTITIES):
        for col, sweep in enumerate(sweeps):
            ax = axes[row, col]
            ax.axhspan(-1.0, 1.0, color=BAND, zorder=0, lw=0)
            ax.axhline(0.0, color=MUTED, lw=0.6, zorder=1)
            for a in atmospheres:
                mine = sorted(
                    (p for p in points if p.sweep == sweep and p.atmosphere == a),
                    key=lambda p: p.value,
                )
                if not mine:
                    continue
                _, colour, marker = STYLE[a]
                x = np.array([p.value for p in mine])
                y = np.array([p.relative(q) * 100 for p in mine])
                shown = np.clip(y, -LIMIT, LIMIT)
                ax.plot(x, shown, color=colour, lw=1.0, zorder=2)
                for xi, yi, si, p in zip(x, y, shown, mine, strict=True):
                    clipped = yi != si
                    ax.plot(
                        xi,
                        si,
                        linestyle="",
                        marker=("^" if yi > 0 else "v") if clipped else marker,
                        markersize=4.5,
                        color=colour,
                        alpha=1.0 if _meaningful(p, q) else 0.3,
                        zorder=3,
                    )
            ax.set_ylim(-LIMIT - 0.3, LIMIT + 0.3)
            ax.set_yticks([-3, -2, -1, 0, 1, 2, 3])
            if row == 0:
                ax.set_title(SWEEP_LABELS[sweep], fontsize=9)
            if col == 0:
                ax.set_ylabel(f"{q} (%)")
    fig.tight_layout(rect=(0, 0.02, 1, 1))
    fig.text(
        0.5,
        0.0,
        "(adjeff − SOS) / SOS.  Shaded: ±1 %.  Faded: SOS below "
        f"{FLOOR}, where its printed decimals limit the comparison.  "
        f"Triangles: beyond ±{LIMIT:g} %, drawn on the edge.",
        ha="center",
        va="bottom",
        color=MUTED,
        fontsize=8,
    )
    _legend(fig, atmospheres, filled=True)
    return fig


@_styled
def distribution(points: list[Point]) -> Figure:
    """Return the relative differences of each quantity, over every point.

    Only points where SOS gives more than `FLOOR` are kept.
    """
    atmospheres = [a for a in STYLE if any(p.atmosphere == a for p in points)]
    fig, ax = plt.subplots(figsize=(8, 3.6))
    ax.axhspan(-1.0, 1.0, color=BAND, zorder=0, lw=0)
    ax.axhline(0.0, color=MUTED, lw=0.6, zorder=1)
    rng = np.random.default_rng(0)
    for i, q in enumerate(QUANTITIES):
        kept = [p for p in points if _meaningful(p, q)]
        ax.boxplot(
            [p.relative(q) * 100 for p in kept],
            positions=[i],
            widths=0.5,
            showfliers=False,
            boxprops={"color": MUTED},
            medianprops={"color": INK, "lw": 1.3},
            whiskerprops={"color": MUTED},
            capprops={"color": MUTED},
            zorder=2,
        )
        for a in atmospheres:
            mine = [p.relative(q) * 100 for p in kept if p.atmosphere == a]
            ax.plot(
                i + rng.uniform(-0.18, 0.18, len(mine)),
                mine,
                zorder=3,
                **_hollow(a, 4.0),
            )
    ax.set_xticks(range(len(QUANTITIES)), QUANTITIES)
    ax.tick_params(axis="x", which="minor", bottom=False, top=False)
    ax.set_ylabel("(adjeff − SOS) / SOS (%)")
    fig.tight_layout()
    _legend(fig, atmospheres)
    return fig


def main() -> None:
    """Draw the three figures of a ``report.json`` into a directory."""
    parser = argparse.ArgumentParser(description="Figures of the SOS comparison.")
    parser.add_argument("--report", type=Path, default=Path("results/report.json"))
    parser.add_argument("--out", type=Path, default=Path("results/figures"))
    args = parser.parse_args()

    points = load(args.report)
    args.out.mkdir(parents=True, exist_ok=True)
    for name, draw in (
        ("sos_vs_adjeff", scatter),
        ("relative_by_sweep", by_sweep),
        ("relative_distribution", distribution),
    ):
        path = args.out / f"{name}.png"
        draw(points).savefig(path, dpi=200, bbox_inches="tight")
        print(f"wrote {path}")


if __name__ == "__main__":
    main()
