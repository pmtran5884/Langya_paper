#!/usr/bin/env python
"""Figure 1C: the volcano, straight from Apron's own gene table.

This script re-scores nothing. Every log-fold change and -log10(P) plotted here
is the portal's. It only

  * reads the direction Apron selected, from whether the mean of its ascending
    (depletion) or descending (enrichment) -log10(p) column is larger,
  * converts the reported average -log10(p) to a p-value,
  * adds a Benjamini-Hochberg FDR within each direction - an ADDITION, since
    the GPP hypergeometric output carries no multiple-testing correction,
  * flags the non-targeting controls, submitted under the portal's NO_SITE_
    naming convention,
  * writes the gene table and renders the figure.

Default canvas is the paper spec: 210 x 130 mm, Arial, every text element at
10 pt. Gene labels are picked separately on each side of zero fold change
(--n-label-pos / --n-label-neg) by Apron score, not as a single top-N.

Usage:
  scripts/04_volcano_apron.py                      # uses apron_submission/
  scripts/04_volcano_apron.py --apron-out apron_out --decorated
"""
import argparse
import glob
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figure_style import apply_figure_style  # noqa: E402

import matplotlib as mpl  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.transforms import Bbox  # noqa: E402

PRETTY = {"LayG1st": "LayV-G round 1", "LayG2nd": "LayV-G round 2"}
RED = "#c8102e"
MM = 25.4


def mean_of_list(s):
    return np.array([np.mean([float(x) for x in str(v).split(";") if x not in ("", "nan")])
                     if isinstance(v, str) else np.nan for v in s])


def bh(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    q = np.empty_like(p)
    m = len(p)
    prev = 1.0
    for rank, idx in enumerate(o[::-1]):
        prev = min(prev, p[idx] * m / (m - rank))
        q[idx] = prev
    return q


def load(path):
    a = pd.read_csv(path, sep="\t")
    a = a.rename(columns={"Gene Symbol": "gene", "Average LFC": "mean_lfc",
                          "Average -log(p-values)": "apron_score",
                          "Number of perturbations": "n_perturbations"})
    asc = mean_of_list(a["Individual ascending -log(p-values)"])
    desc = mean_of_list(a["Individual descending -log(p-values)"])
    a["direction"] = np.where(desc >= asc, "enriched", "depleted")
    a["p_value"] = np.power(10.0, -a.apron_score)
    a["is_control"] = a.gene.astype(str).str.startswith(("NO_SITE", "DUMMY"))
    a["fdr_BH_added"] = np.nan
    for d in ("enriched", "depleted"):
        m = a.direction == d
        a.loc[m, "fdr_BH_added"] = bh(a.loc[m, "p_value"].values)
    return a.set_index("gene")


def _spread(items, x_col, y_lo, y_hi):
    """Evenly distribute label anchors over [y_lo, y_hi], ordered by the
    points' own y so leader lines cannot cross. Even spacing (rather than a
    minimum gap) is what keeps labels legibly apart after the figure is
    scaled down for a multi-panel layout."""
    items = sorted(items, key=lambda t: t[2])
    n = len(items)
    if n == 1:
        ys = [0.5 * (y_lo + y_hi)]
    else:
        ys = [y_lo + i * (y_hi - y_lo) / (n - 1) for i in range(n)]
    return [(it, x_col, y) for it, y in zip(items, ys)]


def volcano_panel(ax, g, title, fdr_cut=0.1, n_pos=3, n_neg=3, legend=True,
                  bare=False, dot_scale=1.0, label_span=(0.42, 0.97)):
    """Draw one arm's volcano onto `ax`.

    bare       : no legend, no axis titles, no panel title (tick labels kept,
                 since without them the fold-change scale is unreadable)
    dot_scale  : multiplier on marker AREA
    n_pos/n_neg: genes labelled on the positive and negative fold-change side,
                 chosen by Apron score within each sign
    """
    y, x, ctrl = g.apron_score.values, g.mean_lfc.values, g.is_control.values
    hit = (~ctrl) & (g.direction.values == "enriched") & (g.fdr_BH_added.values < fdr_cut)
    ax.scatter(x[~ctrl & ~hit], y[~ctrl & ~hit], s=7 * dot_scale, c="#c9c9c9",
               linewidths=0)
    ax.scatter(x[ctrl], y[ctrl], s=14 * dot_scale, facecolors="none",
               edgecolors="#8a8a8a", linewidths=0.6, label="non-targeting control")
    ax.scatter(x[hit], y[hit], s=34 * dot_scale, c=RED, linewidths=0,
               label=f"enriched, FDR < {fdr_cut:g}")

    real = g[~g.is_control]
    pos = real[real.mean_lfc > 0].nlargest(n_pos, "apron_score")
    neg = real[real.mean_lfc < 0].nlargest(n_neg, "apron_score")

    ax.margins(0.08)
    x0, x1 = ax.get_xlim()
    ax.set_xlim(x0 - 0.10 * (x1 - x0), x1 + 0.16 * (x1 - x0))
    y0, y1 = ax.get_ylim()
    ax.set_ylim(y0, y1 + ((0.14 if legend else 0.04) * (y1 - y0)))
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    lo, hi = y0 + label_span[0] * (y1 - y0), y0 + label_span[1] * (y1 - y0)

    placed = []
    for sub, side in ((pos, "right"), (neg, "left")):
        items = [(gene, r.mean_lfc, r.apron_score, r.fdr_BH_added)
                 for gene, r in sub.iterrows()]
        if not items:
            continue
        x_col = (x1 - 0.015 * (x1 - x0)) if side == "right" else (x0 + 0.015 * (x1 - x0))
        placed += [(it, xc, yy, side) for it, xc, yy in _spread(items, x_col, lo, hi)]
    for (gene, gx, gy, fdr), xc, yy, side in placed:
        ax.annotate(gene, (gx, gy), xytext=(xc, yy), textcoords="data",
                    va="center", ha="right" if side == "right" else "left",
                    color=RED if fdr < fdr_cut else "#444444",
                    arrowprops=dict(arrowstyle="-", lw=0.4, color="#b0b0b0",
                                    shrinkA=0, shrinkB=3))

    if not bare:
        ax.set_xlabel("mean sgRNA log$_2$ fold change vs Mock")
        ax.set_title(title, loc="left")
    if legend and not bare:
        ax.legend(frameon=False, loc="upper left", handletextpad=0.4,
                  borderpad=0.2, labelspacing=0.4)


def volcano_single(g, title, out, formats, fdr_cut=0.1, n_pos=3, n_neg=3,
                   width_mm=210.0, height_mm=130.0, bare=False, dot_scale=1.0):
    """One arm, one file, at the exact requested physical size."""
    fig, ax = plt.subplots(figsize=(width_mm / MM, height_mm / MM),
                           layout="constrained")
    volcano_panel(ax, g, title, fdr_cut=fdr_cut, n_pos=n_pos, n_neg=n_neg,
                  legend=not bare, bare=bare, dot_scale=dot_scale)
    if not bare:
        ax.set_ylabel("Apron average $-$log$_{10}$ $P$")
    # Save the exact canvas: bbox_inches="tight" (the rcParam default here)
    # would crop the page to the ink and lose the requested physical size.
    box = Bbox([[0, 0], [width_mm / MM, height_mm / MM]])
    for ext in formats:
        fig.savefig(f"{out}.{ext}", dpi=300, bbox_inches=box, pad_inches=0)
    return fig


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--apron-out", default="apron_submission",
                    help="directory holding the portal's gene table(s)")
    ap.add_argument("--arms", nargs="+", default=["LayG1st"])
    ap.add_argument("--resdir", default="results")
    ap.add_argument("--figdir", default="figures")
    ap.add_argument("--fdr-cut", type=float, default=0.1)
    ap.add_argument("--bare", action="store_true", default=True,
                    help="no legend, axis titles or panel title (ticks kept)")
    ap.add_argument("--decorated", dest="bare", action="store_false",
                    help="restore the legend, axis titles and panel title")
    ap.add_argument("--dot-scale", type=float, default=2.0,
                    help="multiplier on marker AREA")
    ap.add_argument("--n-label-pos", type=int, default=3)
    ap.add_argument("--n-label-neg", type=int, default=3)
    ap.add_argument("--formats", nargs="+", default=["png", "pdf", "svg"])
    ap.add_argument("--font", default="Arial")
    ap.add_argument("--font-size", type=float, default=10.0,
                    help="single point size used for every text element")
    ap.add_argument("--width-mm", type=float, default=210.0)
    ap.add_argument("--height-mm", type=float, default=130.0)
    a = ap.parse_args()

    os.makedirs(a.resdir, exist_ok=True)
    os.makedirs(a.figdir, exist_ok=True)
    fs = a.font_size
    apply_figure_style(frame="open", sizes=(fs, fs, fs), font=a.font)
    mpl.rcParams.update({"mathtext.fontset": "custom", "mathtext.rm": a.font,
                         "mathtext.it": f"{a.font}:italic", "mathtext.bf": f"{a.font}:bold",
                         "mathtext.default": "regular",
                         "font.size": fs, "axes.titlesize": fs,
                         "axes.labelsize": fs, "legend.fontsize": fs,
                         "xtick.labelsize": fs, "ytick.labelsize": fs})

    for arm in a.arms:
        hits = sorted(glob.glob(os.path.join(a.apron_out, f"{arm}*.txt")))
        if not hits:
            raise SystemExit(f"[04] no {arm}*.txt in {a.apron_out}")
        t = load(hits[0])
        csv = os.path.join(a.resdir, f"apron_gene_scores_{arm}_vs_Mock.csv")
        t.to_csv(csv)
        n_hit = int(((t.direction == "enriched") & (t.fdr_BH_added < a.fdr_cut)
                     & (~t.is_control)).sum())
        print(f"[04] {arm}: {len(t):,} genes ({int(t.is_control.sum()):,} control) | "
              f"enriched FDR<{a.fdr_cut:g}: {n_hit} | "
              f"PROCR score {t.loc['PROCR', 'apron_score']:.6f} "
              f"p={t.loc['PROCR', 'p_value']:.3g} FDR={t.loc['PROCR', 'fdr_BH_added']:.3g}")
        print(f"[04] gene table -> {csv}")
        out = os.path.join(a.figdir, f"volcano_apron_{arm}")
        volcano_single(t, f"{PRETTY.get(arm, arm)} vs Mock", out, a.formats,
                       fdr_cut=a.fdr_cut, n_pos=a.n_label_pos, n_neg=a.n_label_neg,
                       width_mm=a.width_mm, height_mm=a.height_mm,
                       bare=a.bare, dot_scale=a.dot_scale)
        print(f"[04] figure   -> {out}.{{{','.join(a.formats)}}} "
              f"({a.width_mm:g} x {a.height_mm:g} mm)")


if __name__ == "__main__":
    main()
