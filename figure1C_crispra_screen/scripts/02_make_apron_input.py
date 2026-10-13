#!/usr/bin/env python
"""Turn the sgRNA count matrix into the two files the Broad GPP (Apron) portal takes.

Apron consumes log-fold changes, never counts, so normalisation is done here:

    lognorm(guide) = log2(reads / total_reads_of_that_arm * 1e6 + 1)
    LFC(guide)     = lognorm(treatment) - lognorm(Mock)

Totals are taken over every guide in the library key, before any filtering.

Two deliberate choices, both of which change the gene ranking:

  * Guides with fewer than --min-reference-reads reads in the unselected Mock
    arm are dropped (1,039 guides read exactly zero). Their fold change is
    undefined; Apron has no minimum-reference-reads option, so leaving them in
    ranks library-absent genes above genuinely depleted ones. They stay in the
    .chip file - Apron only scores barcodes present in the data file.
  * Non-targeting guides are renamed NONTARGET_<n> -> NO_SITE_<n>, the portal's
    own convention. Without this Apron treats each as a one-guide gene instead
    of pooling them into no-site control pseudo-genes.

Outputs (to --outdir):
    screen.chip        Barcode Sequence, Gene Symbol            (tab separated)
    screen_input.txt   Construct Barcode + one LFC column per arm
"""
import argparse
import os

import numpy as np
import pandas as pd


def lognorm(counts):
    """log2 CPM + 1, per column, over all guides in the key."""
    return np.log2(counts / counts.sum(axis=0) * 1e6 + 1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--counts", default="data/sgrna_count_matrix.csv",
                    help="CSV with Guide_Seq, Gene and one integer column per arm")
    ap.add_argument("--reference", default="Mock", help="unselected arm")
    ap.add_argument("--treatments", nargs="+", default=["LayG1st", "LayG2nd"])
    ap.add_argument("--min-reference-reads", type=int, default=1)
    ap.add_argument("--outdir", default="apron_submission")
    ap.add_argument("--chip-name", default="layg_crispra.chip")
    ap.add_argument("--data-name", default="layg_crispra_input.txt")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    m = pd.read_csv(a.counts)
    arms = [a.reference] + a.treatments
    counts = m[arms].astype(np.int64)

    # --- .chip : every barcode in the library, controls under NO_SITE_
    sym = m["Gene"].astype(str).str.replace(r"^NONTARGET_", "NO_SITE_", regex=True)
    chip = pd.DataFrame({"Barcode Sequence": m["Guide_Seq"], "Gene Symbol": sym})
    chip_path = os.path.join(a.outdir, a.chip_name)
    chip.to_csv(chip_path, sep="\t", index=False)

    # --- data file : one log-fold-change column per treatment arm
    ln = lognorm(counts)
    data = pd.DataFrame({"Construct Barcode": m["Guide_Seq"]})
    for t in a.treatments:
        data[f"{t}_vs_{a.reference}"] = ln[t] - ln[a.reference]
    keep = counts[a.reference] >= a.min_reference_reads
    data = data[keep]
    data_path = os.path.join(a.outdir, a.data_name)
    data.to_csv(data_path, sep="\t", index=False, float_format="%.9f")

    n_ctrl = int(sym.str.startswith("NO_SITE").sum())
    avg = chip[~sym.str.startswith("NO_SITE")].groupby("Gene Symbol").size().mean()
    print(f"[02] {chip_path}: {len(chip):,} barcodes, {n_ctrl:,} no-site controls, "
          f"{avg:.2f} guides/gene")
    print(f"[02] {data_path}: {len(data):,} barcodes x {len(a.treatments)} arm(s) "
          f"({int((~keep).sum()):,} dropped for < {a.min_reference_reads} "
          f"{a.reference} read(s))")
    print("[02] portal settings used for the paper: Hypergeometric, directionality "
          "Both, top 100%, dummy gene set size 4, min/max guides 2/4, "
          "no-site controls displayed")


if __name__ == "__main__":
    main()
