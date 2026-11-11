"""Figures of the note's section on realizations, from the results of fca.scans.run, written to
fca/note/figures. All of them, or those named, from the repository root:

    venv/bin/python -m fca.scans.figures [lifetime angle aspect tilt spread spacing ensemble]
"""

import json
import sys
from dataclasses import replace
from pathlib import Path

import figure_provenance
import matplotlib.pyplot as plt
import numpy as np
from figure_provenance import provenance

from fca.scans.model import Case, expected
from fca.scans.run import RESULTS, SPREAD_LIFETIME

plt.style.use(["cosmoplots.default"])
plt.rcParams["text.latex.preamble"] = r"\usepackage{amsmath} \usepackage{mathptmx} \usepackage{bm}"
FONTSIZE = 11  # the note's 11 pt preprint text, in Times as mathptmx sets it
plt.rcParams.update({"font.family": "serif", "font.serif": ["Times"], "font.size": FONTSIZE,
                     "axes.labelsize": FONTSIZE, "legend.fontsize": 9, "xtick.labelsize": 9, "ytick.labelsize": 9})

TEXTWIDTH = 468 / 72.27  # inches, revtex4-2 preprint on a4
OUT = Path(__file__).parents[1] / "note" / "figures"
GREY = "0.55"

# Okabe-Ito colour, marker, line style of the model line, and label of each method
STYLE = {
    "2TDE": ("#CC79A7", "v", "-.", "2TDE"),
    "3TDE": ("#E69F00", "o", "-", "3TDE"),
    "CROSSING": ("#0072B2", "D", "--", "crossing time"),
    "FCA": ("#009E73", "s", ":", "FCA"),
}
# quantity, methods, label, value without error
PANELS = (
    ("x", ("2TDE", "CROSSING"), r"$\widehat{v}_x/v_x$, two points", 1),
    ("speed", ("3TDE", "CROSSING", "FCA"), r"$|\widehat{\bm{v}}|/|\bm{v}|$, three points", 1),
    ("direction", ("3TDE", "CROSSING", "FCA"), r"direction of $\widehat{\bm{v}}$ from $\bm{v}$ (degrees)", 0),
)


def load(name):
    """[(Case, {method: array (2,)})], one item per realization of the scan."""
    rows = json.loads((RESULTS / f"{name}.json").read_text())["rows"]
    fields = Case.__dataclass_fields__
    return [(Case(**{k: (np.inf if r[k] is None else r[k]) for k in fields}),
             {m: np.array(v, dtype=float) for m, v in r["estimates"].items()}) for r in rows]


def errors(case, estimates):
    """{quantity: {method: value}} of one set of estimates: the x component over the true one, the
    speed over the true one, and the direction from the true one in degrees."""
    v = case.velocity()
    out = {"x": {}, "speed": {}, "direction": {}, "vector": {}, "component": {}}
    for method, est in estimates.items():
        out["x"][method] = est[0] / v[0]
        out["speed"][method] = np.hypot(*est)
        turn = np.degrees(np.arctan2(est[1], est[0]) - np.arctan2(v[1], v[0]))
        out["direction"][method] = (turn + 180) % 360 - 180
        out["vector"][method] = np.hypot(*(est - v))
        out["component"][method] = abs(est[0] - v[0])
    return out


def draw_scan(axes, results, base, parameter, grid, xscale="linear", infinite_at=None):
    """One row of PANELS against a parameter of the case. Markers: mean and standard deviation over
    the records of each case of `results` that equals `base` but for the parameter. Lines: what
    the method gives on the model correlation function. Grey: no error. The methods are moved
    apart along the abscissa so that their bars do not overlap. A parameter value of infinity is
    drawn at infinite_at."""
    values = sorted({getattr(c, parameter) for c, _ in results if replace(c, **{parameter: getattr(base, parameter)}) == base})
    model = [errors(c, expected(c)) for c in (replace(base, **{parameter: x}) for x in grid)]
    place = lambda x: infinite_at if np.isinf(x) else x
    for ax, (quantity, methods, label, reference) in zip(axes, PANELS):
        ax.axhline(reference, color=GREY, lw=0.8, zorder=1)
        for j, method in enumerate(methods):
            color, marker, ls, _ = STYLE[method]
            ax.plot(grid, [m[quantity][method] for m in model], color=color, ls=ls, lw=1.0, zorder=2)
            shift = (j - (len(methods) - 1) / 2) * 0.016
            for x in values:
                case = replace(base, **{parameter: x})
                reads = np.array([errors(c, e)[quantity][method] for c, e in results if c == case])
                reads = reads[np.isfinite(reads)]
                if len(reads) < 2:
                    continue
                span = np.ptp(ax.get_xlim()) if xscale == "linear" else None
                xs = place(x) + shift * span if span else place(x) * 10 ** (shift * 1.6)
                ax.errorbar(xs, reads.mean(), reads.std(ddof=1), fmt=marker, color=color, ms=3.5, mew=0.8,
                            mfc="w" if method in ("2TDE", "CROSSING") else color, elinewidth=0.8, capsize=0, zorder=4)
        ax.set_ylabel(label)
        ax.set_xscale(xscale)


def finish(fig):
    """Tighten the layout and put one legend of the four methods above the panels, where it covers
    no data."""
    fig.tight_layout(w_pad=0.6)
    handles = [plt.Line2D([], [], color=STYLE[m][0], marker=STYLE[m][1], ls=STYLE[m][2], lw=1.0, ms=3.5, mew=0.8,
                          mfc="w" if m in ("2TDE", "CROSSING") else STYLE[m][0]) for m in STYLE]
    fig.legend(handles, [STYLE[m][3] for m in STYLE], handlelength=2.6, ncol=4, loc="lower center",
               bbox_to_anchor=(0.5, 1.0), frameon=False)
    return fig


def row(n_rows=1, height=2.35):
    fig, axes = plt.subplots(n_rows, 3, figsize=(TEXTWIDTH, height * n_rows), squeeze=False)
    for ax, letter in zip(axes.flat, "abcdef"):
        ax.set_title(f"({letter})", loc="left", fontsize=FONTSIZE)
    return fig, axes


@provenance
def fig_lifetime():
    """Decorrelation: round pulses moving along x, against their lifetime (the frozen pattern at the
    right end). (a) The two-point estimates of the pair along x, (b) the speed and (c) the direction
    of the three-point estimates, as ratios to and angles from the true velocity. Markers: mean
    and standard deviation over eight records; lines: the model. Logarithmic ordinates in (a)
    and (b), where the time delays overestimate by up to a factor of five."""
    fig, axes = row()
    for ax in axes[0]:
        ax.set_xlim(0.35, 60)
    draw_scan(axes[0], load("lifetime"), Case(), "taup", np.geomspace(0.4, 22, 60), "log", infinite_at=40)
    for ax in axes[0]:
        ax.set_xlabel(r"$\tau_\text{p} v/\ell$")
        ax.set_xticks([1, 10, 40], ["1", "10", r"$\infty$"])
    for ax in axes[0, :2]:
        ax.set_yscale("log")
        ax.set_ylim(0.7, 7)
        ax.set_yticks([1, 2, 5], ["1", "2", "5"])
        ax.minorticks_off()
    axes[0, 2].set_ylim(-25, 25)
    return finish(fig)


@provenance
def fig_angle():
    """Misalignment: round frozen pulses, against the angle of their velocity from x. Panels and
    symbols as in fig_lifetime. The two-point delay along x gives v_x / cos^2(beta); the other
    estimates have no bias in the model. The crossing time and the FCA scatter and fall short in
    the component across the motion where one pair has no crossing (beta = 0) or a late one
    (beta = 15 and 75 degrees), in the tails of the correlation functions."""
    fig, axes = row()
    for ax in axes[0]:
        ax.set_xlim(-6, 81)
    draw_scan(axes[0], load("angle"), Case(), "angle", np.linspace(0, 78, 53))
    for ax in axes[0]:
        ax.set_xlabel(r"$\beta$ (degrees)")
        ax.set_xticks([0, 30, 60])
    axes[0, 0].set_yscale("log")
    axes[0, 0].set_ylim(0.7, 25)
    axes[0, 0].set_yticks([1, 2, 5, 10, 20], ["1", "2", "5", "10", "20"])
    axes[0, 0].minorticks_off()
    axes[0, 1].set_ylim(0.8, 1.2)
    axes[0, 2].set_ylim(-25, 25)
    return finish(fig)


@provenance
def fig_aspect():
    """Barber pole against the elongation: frozen pulses moving along x with their long axis 30
    degrees from x, against the aspect ratio. Panels and symbols as in fig_lifetime. The pair
    along x has no bias for a frozen pattern moving along it; the 3TDE and the crossing time turn
    towards the short axis."""
    fig, axes = row()
    for ax in axes[0]:
        ax.set_xlim(0.8, 4.2)
    draw_scan(axes[0], load("aspect"), Case(tilt=30), "aspect", np.linspace(1, 4.1, 63))
    for ax in axes[0]:
        ax.set_xlabel(r"$\ell_1/\ell_2$")
        ax.set_xticks([1, 2, 3, 4])
    axes[0, 0].set_ylim(0.8, 1.2)
    axes[0, 1].set_ylim(0, 1.6)
    axes[0, 2].set_ylim(-75, 25)
    return finish(fig)


@provenance
def fig_tilt(aspects=(2.0, 4.0)):
    """Barber pole against the tilt: frozen pulses moving along x, against the angle of their long
    axis from x, at the aspect ratios 2 (top row) and 4 (bottom row). Panels and symbols as in
    fig_lifetime. One row per aspect ratio so that each panel holds the methods only."""
    fig, axes = row(2)
    results = load("tilt")
    for axes_row, aspect in zip(axes, aspects):
        for ax in axes_row:
            ax.set_xlim(-8, 98)
        draw_scan(axes_row, results, Case(aspect=aspect), "tilt", np.linspace(0, 90, 91))
        for ax in axes_row:
            ax.set_xticks([0, 30, 60, 90])
        axes_row[0].set_ylim(0.8, 1.2)
        axes_row[1].set_ylim(0, 2.4)
        axes_row[2].set_ylim(-90, 30)
        axes_row[2].text(0.96, 0.06, rf"$\ell_1/\ell_2 = {aspect:g}$", transform=axes_row[2].transAxes, ha="right", fontsize=9)
    for ax in axes[1]:
        ax.set_xlabel(r"$\alpha$ (degrees)")
    return finish(fig)


@provenance
def fig_spread():
    """A spread of velocities, which none of the estimators models: pulses with the lifetime
    SPREAD_LIFETIME moving along x on average, against the standard deviation of each velocity
    component, for round pulses (top row) and for pulses of aspect ratio 2 tilted by 30 degrees
    (bottom row). Panels and symbols as in fig_lifetime; the lines are the estimators on the model
    correlation with the lag-dependent shape W + sigma^2 tau^2."""
    fig, axes = row(2)
    for axes_row, name, base, text in zip(
        axes, ("spread_round", "spread_tilted"),
        (Case(taup=SPREAD_LIFETIME), Case(taup=SPREAD_LIFETIME, aspect=2, tilt=30)),
        ("round", r"$\ell_1/\ell_2 = 2$, $\alpha = 30^\circ$"),
    ):
        for ax in axes_row:
            ax.set_xlim(-0.05, 0.55)
        draw_scan(axes_row, load(name), base, "spread", np.linspace(0, 0.52, 53))
        for ax in axes_row:
            ax.set_xticks([0, 0.2, 0.4])
        axes_row[0].set_ylim(0.6, 2.0)
        axes_row[1].set_ylim(0.6, 2.0)
        axes_row[2].set_ylim(-60, 30)
        axes_row[2].text(0.96, 0.06, text, transform=axes_row[2].transAxes, ha="right", fontsize=9)
    for ax in axes[1]:
        ax.set_xlabel(r"$\sigma_v/v$")
    return finish(fig)


@provenance
def fig_spacing():
    """The pixel separation: pulses of aspect ratio 2 tilted by 30 degrees, with the lifetime 2,
    moving along x, against the separation over the pulse size. Panels and symbols as in
    fig_lifetime. The biases of the time delays and of the crossing time do not depend on the
    separation; the scatter of the FCA grows at small separations, and at the largest one the
    shape fit fails (no marker)."""
    fig, axes = row()
    for ax in axes[0]:
        ax.set_xlim(0.28, 2.5)
    draw_scan(axes[0], load("spacing"), Case(taup=2, aspect=2, tilt=30), "spacing", np.geomspace(0.3, 2.3, 50), "log")
    for ax in axes[0]:
        ax.set_xlabel(r"$\Delta/\ell$")
        ax.set_xticks([0.5, 1, 2], ["0.5", "1", "2"])
        ax.set_xticks([], minor=True)
    axes[0, 0].set_ylim(0.6, 2.0)
    axes[0, 1].set_ylim(0, 2.0)
    axes[0, 2].set_ylim(-75, 50)
    return finish(fig)


@provenance
def fig_ensemble():
    """Every parameter drawn at random (fca.scans.run.ensemble_cases), one record per case: the
    fraction of cases with an error below the abscissa. (a) The component along x from two points,
    |v_x estimate - v_x| / |v|; (b) the vector from three points, |v estimate - v| / |v|. Thick
    lines with markers: the realizations; thin lines: the estimators on the model correlation
    function, the bias alone. Cumulative distributions on a logarithmic abscissa so that errors
    from a hundredth of the speed to several times it show together."""
    results = load("ensemble")
    fig, axes = plt.subplots(1, 2, figsize=(TEXTWIDTH, 2.7))
    for ax, letter, quantity, methods, label in zip(
        axes, "ab", ("component", "vector"), (PANELS[0][1], PANELS[1][1]),
        (r"$|\widehat{v}_x - v_x|/|\bm{v}|$, two points", r"$|\widehat{\bm{v}} - \bm{v}|/|\bm{v}|$, three points"),
    ):
        for method in methods:
            color, marker, ls, name = STYLE[method]
            reads = np.array([errors(c, e)[quantity][method] for c, e in results])
            model = np.array([errors(c, expected(c))[quantity][method] for c, _ in results])
            fraction = np.arange(1, len(reads) + 1) / len(reads)
            # a case without an estimate counts as an error beyond the axis
            ax.plot(np.sort(np.where(np.isfinite(model), model, np.inf)), fraction, color=color, ls=ls, lw=0.7)
            ax.plot(np.sort(np.where(np.isfinite(reads), reads, np.inf)), fraction, color=color, ls=ls, lw=1.5,
                    marker=marker, ms=3.5, mew=0.8, mfc="w" if method in ("2TDE", "CROSSING") else color,
                    markevery=(0.05, 0.12), label=name)
        ax.set_xscale("log")
        ax.set_xlim(3e-3, 10)
        ax.set_ylim(0, 1)
        ax.set_xlabel(label)
        ax.set_title(f"({letter})", loc="left", fontsize=FONTSIZE)
        ax.legend(loc="upper left", handlelength=2.6)
    axes[0].set_ylabel("fraction of cases")
    fig.tight_layout(w_pad=1.0)
    return fig


def save(fig, name):
    figure_provenance.save(fig, OUT / f"{name}.pdf", stacklevel=2, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)


FIGURES = {"lifetime": fig_lifetime, "angle": fig_angle, "aspect": fig_aspect, "tilt": fig_tilt,
           "spread": fig_spread, "spacing": fig_spacing, "ensemble": fig_ensemble}


def main(names):
    OUT.mkdir(exist_ok=True)
    for name in names or FIGURES:
        save(FIGURES[name](), f"scan_{name}")


if __name__ == "__main__":
    main(sys.argv[1:])
