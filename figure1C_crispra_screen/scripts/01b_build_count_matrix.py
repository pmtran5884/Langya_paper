#!/usr/bin/env python
"""Join the per-sample count files from step 01 into one count matrix.

Step 01 writes `<sample>_counts.tsv` as a sparse (guide_idx, count) pair list,
indexed by row order in the library key. This expands those back onto the full
key - every guide present, zeros included - and writes the matrix that steps 02
and onward read.

Usage:
  scripts/01b_build_count_matrix.py --key data/formatted_sgrna_key.csv \
      --countdir counts --samples Mock LayG1st LayG2nd \
      --out data/sgrna_count_matrix.csv
"""
import argparse
import os

import numpy as np
import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--key", default="data/formatted_sgrna_key.csv")
    ap.add_argument("--countdir", default="counts")
    ap.add_argument("--samples", nargs="+", default=["Mock", "LayG1st", "LayG2nd"])
    ap.add_argument("--out", default="data/sgrna_count_matrix.csv")
    a = ap.parse_args()

    key = pd.read_csv(a.key)
    mat = pd.DataFrame(index=range(len(key)))
    for s in a.samples:
        path = os.path.join(a.countdir, f"{s}_counts.tsv")
        c = pd.read_csv(path, sep="\t")
        v = np.zeros(len(key), dtype=np.int64)
        v[c.guide_idx.values] = c["count"].values
        mat[s] = v
        print(f"[01b] {s}: {v.sum():,} mapped reads over {(v > 0).sum():,} guides")
    mat["Gene"] = key["Gene"].values
    mat["Guide_Seq"] = key["Guide_Seq"].values
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    mat.to_csv(a.out, index=False)
    print(f"[01b] {a.out}: {len(mat):,} guides x {len(a.samples)} arms")


if __name__ == "__main__":
    main()
