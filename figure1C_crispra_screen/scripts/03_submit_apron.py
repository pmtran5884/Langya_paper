#!/usr/bin/env python
"""Submit the screen to the Broad GPP pooled-screen portal (Apron) and collect it.

The portal has no documented API; this drives the public web form with a
multipart POST and then polls the results page for the download links. The form
fields below are the ones the paper used:

    analysis type        Hypergeometric   (exact string; a lowercase value
                                           silently falls through to STARS)
    directionality       both             (so depletion is in the same table)
    % of guides included 100              (score all of a gene's guides)
    dummy gene set size  4                (matches 3.97 mean guides/gene)
    guides per gene      min 2, max 4     (excludes genes left with one guide)
    no-site controls     displayed

Jobs finish in ~2-4 min. The portal returns one gene table (.txt) per LFC
column in the data file, plus its own volcano PDFs.

Usage:
  scripts/03_submit_apron.py --chip apron_submission/layg_crispra.chip \
      --data apron_submission/layg_crispra_input.txt --outdir apron_out

Nothing downstream depends on re-running this: apron_submission/ already holds
the exact files submitted for the paper and the portal's unmodified output, so
step 04 runs offline. A reviewer can also upload the two files by hand at
https://portals.broadinstitute.org/gpp/public/analysis-tools/crispr-gene-scoring
with the settings above and get the same table.
"""
import argparse
import html
import os
import re
import time

import requests

URL = "https://portals.broadinstitute.org/gpp/public/analysis-tools/crispr-gene-scoring"
RESULTS = URL + "-results"
BASE = "https://portals.broadinstitute.org"


def chip_meta(chip_path):
    """Barcode count, no-site count and mean guides/gene, as the form wants them."""
    n, nosite, per_gene = 0, 0, {}
    with open(chip_path) as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 2:
                continue
            n += 1
            if f[1].startswith("NO_SITE"):
                nosite += 1
            else:
                per_gene[f[1]] = per_gene.get(f[1], 0) + 1
    avg = sum(per_gene.values()) / len(per_gene) if per_gene else 0.0
    return dict(n_barcodes=n, n_nosite=nosite, avg_guides=avg)


def submit(chip_path, data_path, meta, directionality="both", top_pct=100,
           dummy_size=4, min_guides=2, max_guides=4, genes_to_label=5):
    s = requests.Session()
    s.headers.update({"Referer": URL})
    files = {"chipfile-input": (os.path.basename(chip_path), open(chip_path, "rb"), "text/plain"),
             "datafile-input": (os.path.basename(data_path), open(data_path, "rb"), "text/plain")}
    form = {"clonepool_type": "strict", "clonepool_json_strict": "", "clonepool_json_lax": "",
            "analysistype": "Hypergeometric", "directionality": directionality,
            "percentincluded": str(top_pct), "dummygenesetsize": str(dummy_size),
            "barcodespergenemin": str(min_guides), "barcodespergenemax": str(max_guides),
            "genestohighlight": str(genes_to_label), "displaynositecontrols": "on",
            "numbarcodes": str(meta["n_barcodes"]),
            "numbarcodespergene": f"{meta['avg_guides']:.2f}",
            "numnositebarcodes": str(meta["n_nosite"]),
            "numoneintergenicsitebarcodes": "0",
            "threshold": "10", "includefirst": "Yes", "userclickedsubmit": "true"}
    r = s.post(URL, data=form, files=files, timeout=600)
    m = re.search(r"reqid=([^&\"\s]+)", r.url)
    if not m:
        raise SystemExit("[03] the portal did not return a job id")
    return s, m.group(1).replace("%3A", ":")


def poll(session, reqid, outdir, timeout_s=5400, interval_s=30):
    os.makedirs(outdir, exist_ok=True)
    deadline = time.time() + timeout_s
    while time.time() < deadline:
        r = session.get(RESULTS, params={"reqid": reqid}, timeout=120)
        # Exclude only the portal's own documentation/example downloads, by
        # their actual names. A bare "sample" substring test would also throw
        # away real results whenever an arm is named e.g. "sample_1".
        links = [l for l in sorted(set(re.findall(r'href="([^"]*download[^"]*)"', r.text)))
                 if "README" not in l and "crispr_gene_scoring_sample" not in l]
        if links:
            saved = []
            for l in links:
                u = l if l.startswith("http") else BASE + l
                fn = re.search(r"filename=([^&]+)", u)
                fn = html.unescape(fn.group(1)) if fn else os.path.basename(u.split("?")[0])
                g = session.get(u, timeout=600)
                if g.ok and len(g.content) > 100:
                    path = os.path.join(outdir, fn)
                    with open(path, "wb") as fh:
                        fh.write(g.content)
                    saved.append(path)
                    print(f"[03] saved {path} ({len(g.content):,} bytes)")
            if saved:
                return saved
        print(f"[03] job {reqid} processing ...", flush=True)
        time.sleep(interval_s)
    raise SystemExit("[03] timed out waiting for the portal")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chip", default="apron_submission/layg_crispra.chip")
    ap.add_argument("--data", default="apron_submission/layg_crispra_input.txt")
    ap.add_argument("--outdir", default="apron_out")
    ap.add_argument("--reqid", default=None,
                    help="skip submission and poll an existing job id")
    a = ap.parse_args()
    if a.reqid:
        s = requests.Session()
        s.headers.update({"Referer": URL})
        poll(s, a.reqid, a.outdir)
        return
    meta = chip_meta(a.chip)
    print(f"[03] submitting {meta['n_barcodes']:,} barcodes "
          f"({meta['n_nosite']:,} no-site, {meta['avg_guides']:.2f} guides/gene)")
    s, reqid = submit(a.chip, a.data, meta)
    print(f"[03] job id {reqid}")
    poll(s, reqid, a.outdir)


if __name__ == "__main__":
    main()
