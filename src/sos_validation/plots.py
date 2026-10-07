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
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

from .five_s import QUANTITIES  # noqa: E402
from .report import FLOOR, Point  # noqa: E402

#: Atmosphere to ``(label, colour, marker)``.  Colour and shape both
#: carry identity, so the figures read in greyscale and for any vision.
STYLE = {
    "rayleigh": ("Molecules only", "#2a78d6", "o"),
    "aerosol": ("Aerosol only", "#eb6834", "s"),
    "both": ("Molecules and aerosol", "#1baf7a", "D"),
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
GRID = "#e4e3df"
BAND = "#efeeea"

plt.rcParams.update(
    {
        "figure.facecolor": "#fcfcfb",
        "axes.facecolor": "#fcfcfb",
        "axes.edgecolor": MUTED,
        "axes.labelcolor": INK,
        "axes.titlecolor": INK,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.color": GRID,
        "grid.linewidth": 0.8,
        "xtick.color": MUTED,
        "ytick.color": MUTED,
        "font.size": 10,
        "savefig.dpi": 150,
        "savefig.bbox": "tight",
    }
)


def load(path: Path) -> list[Point]:
    """Return the points of a ``report.json``."""
    return [Point(**entry) for entry in json.loads(path.read_text())]


def _meaningful(point: Point, quantity: str) -> bool:
    """Return whether SOS is far enough from zero for a relative difference."""
    return abs(point.sos[quantity]) > FLOOR


def _legend(fig: Figure, atmospheres: list[str]) -> None:
    """Add the atmosphere legend above the panels."""
    handles = [
        Line2D(
            [],
            [],
            linestyle="",
            marker=STYLE[a][2],
            color=STYLE[a][1],
            markeredgecolor="white",
            markersize=8,
            label=STYLE[a][0],
        )
        for a in atmospheres
    ]
    fig.legend(
        handles=handles,
        loc="upper center",
        ncol=len(handles),
        frameon=False,
        bbox_to_anchor=(0.5, 1.0),
    )


def scatter(points: list[Point]) -> Figure:
    """Return adjeff against SOS, one panel per quantity, with ``y = x``."""
    atmospheres = [a for a in STYLE if any(p.atmosphere == a for p in points)]
    fig, axes = plt.subplots(2, 3, figsize=(12, 8.2))
    for ax, q in zip(axes.flat, QUANTITIES, strict=True):
        sos = np.array([p.sos[q] for p in points])
        adj = np.array([p.adjeff[q] for p in points])
        low, high = min(sos.min(), adj.min()), max(sos.max(), adj.max())
        pad = 0.05 * (high - low)
        ax.plot(
            [low - pad, high + pad], [low - pad, high + pad], color=MUTED, lw=1, ls="--"
        )
        for a in atmospheres:
            mine = [p for p in points if p.atmosphere == a]
            ax.scatter(
                [p.sos[q] for p in mine],
                [p.adjeff[q] for p in mine],
                marker=STYLE[a][2],
                s=42,
                color=STYLE[a][1],
                edgecolors="white",
                linewidths=1.0,
                zorder=3,
            )
        ax.set_xlim(low - pad, high + pad)
        ax.set_ylim(low - pad, high + pad)
        ax.set_aspect("equal")
        ax.set_title(f"`{q}`".strip("`"), fontsize=11, loc="left")
        ax.set_xlabel("SOS")
        ax.set_ylabel("adjeff")

        rel = [abs(p.relative(q)) for p in points if _meaningful(p, q)]
        if rel:
            ax.text(
                0.97,
                0.04,
                f"median {np.median(rel) * 100:.2f} %\nworst {max(rel) * 100:.2f} %",
                transform=ax.transAxes,
                ha="right",
                va="bottom",
                color=MUTED,
                fontsize=9,
            )
    _legend(fig, atmospheres)
    fig.suptitle(
        "adjeff against SOS on the six 5S terms",
        y=1.04,
        fontsize=13,
        color=INK,
    )
    fig.tight_layout()
    return fig


def by_sweep(points: list[Point]) -> Figure:
    """Return the relative difference against each swept parameter.

    Rows are quantities, columns sweeps.  Hollow markers are points where
    SOS gives less than `FLOOR`, whose relative difference means little.
    """
    sweeps = [s for s in SWEEP_LABELS if any(p.sweep == s for p in points)]
    atmospheres = [a for a in STYLE if any(p.atmosphere == a for p in points)]
    fig, axes = plt.subplots(
        len(QUANTITIES),
        len(sweeps),
        figsize=(2.9 * len(sweeps), 1.9 * len(QUANTITIES)),
        sharey=True,
        sharex="col",
        squeeze=False,
    )
    for row, q in enumerate(QUANTITIES):
        for col, sweep in enumerate(sweeps):
            ax = axes[row, col]
            ax.axhspan(-1.0, 1.0, color=BAND, zorder=0, lw=0)
            ax.axhline(0.0, color=MUTED, lw=0.8, zorder=1)
            for a in atmospheres:
                mine = sorted(
                    (p for p in points if p.sweep == sweep and p.atmosphere == a),
                    key=lambda p: p.value,
                )
                if not mine:
                    continue
                label, colour, marker = STYLE[a]
                x = np.array([p.value for p in mine])
                y = np.array([p.relative(q) * 100 for p in mine])
                shown = np.clip(y, -LIMIT, LIMIT)
                ax.plot(x, shown, color=colour, lw=1.5, zorder=2)
                for xi, yi, si, p in zip(x, y, shown, mine, strict=True):
                    clipped = yi != si
                    ax.plot(
                        xi,
                        si,
                        linestyle="",
                        marker=("^" if yi > 0 else "v") if clipped else marker,
                        markersize=6.5,
                        color=colour,
                        markerfacecolor=colour if _meaningful(p, q) else "white",
                        markeredgecolor=colour,
                        zorder=3,
                    )
            ax.set_ylim(-LIMIT - 0.3, LIMIT + 0.3)
            ax.set_yticks([-3, -2, -1, 0, 1, 2, 3])
            if row == 0:
                ax.set_title(SWEEP_LABELS[sweep], fontsize=10)
            if col == 0:
                ax.set_ylabel(f"{q}\n(%)")
    _legend(fig, atmospheres)
    fig.suptitle(
        "Relative difference (adjeff − SOS) / SOS, against each swept parameter",
        y=1.02,
        fontsize=13,
        color=INK,
    )
    fig.text(
        0.5,
        -0.01,
        "Shaded: ±1 %.  Hollow: SOS below "
        f"{FLOOR}, where its printed decimals limit the comparison.  "
        f"Triangles: beyond ±{LIMIT:g} %, drawn on the edge.",
        ha="center",
        color=MUTED,
        fontsize=9,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    return fig


def distribution(points: list[Point]) -> Figure:
    """Return the relative differences of each quantity, over every point.

    Only points where SOS gives more than `FLOOR` are kept.
    """
    atmospheres = [a for a in STYLE if any(p.atmosphere == a for p in points)]
    fig, ax = plt.subplots(figsize=(10, 4.6))
    ax.axhspan(-1.0, 1.0, color=BAND, zorder=0, lw=0)
    ax.axhline(0.0, color=MUTED, lw=0.8, zorder=1)
    rng = np.random.default_rng(0)
    every = []
    for i, q in enumerate(QUANTITIES):
        kept = [p for p in points if _meaningful(p, q)]
        values = [p.relative(q) * 100 for p in kept]
        every += [abs(v) for v in values]
        ax.boxplot(
            values,
            positions=[i],
            widths=0.5,
            showfliers=False,
            patch_artist=True,
            boxprops={"facecolor": "white", "edgecolor": MUTED},
            medianprops={"color": INK, "lw": 1.5},
            whiskerprops={"color": MUTED},
            capprops={"color": MUTED},
            zorder=2,
        )
        for a in atmospheres:
            mine = [p.relative(q) * 100 for p in kept if p.atmosphere == a]
            ax.scatter(
                i + rng.uniform(-0.18, 0.18, len(mine)),
                mine,
                marker=STYLE[a][2],
                s=22,
                color=STYLE[a][1],
                edgecolors="white",
                linewidths=0.6,
                zorder=3,
            )
    ax.set_xticks(range(len(QUANTITIES)), QUANTITIES)
    ax.set_ylabel("(adjeff − SOS) / SOS (%)")
    ax.grid(axis="x", visible=False)
    _legend(fig, atmospheres)
    fig.suptitle(
        f"Median |difference| over {len(every)} comparisons: {np.median(every):.2f} %",
        y=1.04,
        fontsize=13,
        color=INK,
    )
    fig.tight_layout()
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
        draw(points).savefig(path)
        print(f"wrote {path}")


if __name__ == "__main__":
    main()
