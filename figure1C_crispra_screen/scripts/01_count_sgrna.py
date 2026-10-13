#!/usr/bin/env python
"""Read-level QC + sgRNA counting for a px52-amplified CRISPRa screen.

Amplicon architecture of R1 (px52 forward primers F1-F4):
    [Illumina adapter] [0-3 nt stagger] CTTGTGGAAAGGACGAAACACCG [20 nt sgRNA] GTTTAAGAGCTATGCTG...
The stagger (F1="", F2="A", F3="TC", F4="ACT") is consumed by the sequencer as
part of read 1, so the sgRNA starts at a variable offset. Reads are therefore
anchored on the constant vector sequence rather than a fixed position.

Outputs per sample (to --outdir):
    <sample>_counts.csv   Guide_Seq, Gene, Guide_Num, count   (exact 20-mer matches)
    <sample>_qc.json      read fates, primer-stagger usage, unmatched 20-mers, R2 entropy
"""
import argparse, collections, gzip, json, os, sys

# constant vector sequence immediately 5' of the sgRNA; ends on the U6 +1 G
ANCHOR = "GGAAAGGACGAAACACCG"
PRIMER_3P = "CTTGTGGAAAGGACGAAACACCG"      # full px52 genome-annealing part
SCAFFOLD = "GTTTAAGAGCTATGCTG"             # tracrRNA, immediately 3' of the sgRNA
STAGGERS = {0: "px52_F1", 1: "px52_F2", 2: "px52_F3", 3: "px52_F4"}
SG_LEN = 20
SCAN_WINDOW = 40        # anchor must start within the first 40 nt of R1


def fastq_seqs(path):
    with gzip.open(path, "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 1:
                yield line.rstrip("\n")


def count_sample(name, r1, r2, key_index, outdir, r2_subsample=200_000):
    fate = collections.Counter()
    stagger = collections.Counter()
    counts = collections.Counter()
    unmatched = collections.Counter()
    lens = collections.Counter()
    n = 0
    for seq in fastq_seqs(r1):
        n += 1
        lens[len(seq)] += 1
        p = seq.find(ANCHOR, 0, SCAN_WINDOW)
        if p < 0:
            # salvage: is this an on-target amplicon with a mutated anchor?
            fate["no_anchor_scaffold_present" if SCAFFOLD[:12] in seq
                 else "no_anchor_no_scaffold"] += 1
            continue
        off = p - (len(PRIMER_3P) - len(ANCHOR))   # inferred stagger length
        stagger[off if 0 <= off <= 3 else -1] += 1
        s = p + len(ANCHOR)
        guide = seq[s:s + SG_LEN]
        if len(guide) < SG_LEN:
            fate["anchor_truncated_read"] += 1
            continue
        scaffold_ok = seq[s + SG_LEN:s + SG_LEN + 10] == SCAFFOLD[:10]
        idx = key_index.get(guide)
        if idx is not None:
            counts[idx] += 1
            fate["mapped_scaffold_ok" if scaffold_ok else "mapped_scaffold_bad"] += 1
        else:
            unmatched[guide] += 1
            fate["unmapped_scaffold_ok" if scaffold_ok else "unmapped_scaffold_bad"] += 1
        if n % 2_000_000 == 0:
            print(f"[{name}] {n/1e6:.0f}M reads", flush=True)

    # R2 positional base composition: is R2 informative (variable) at all?
    r2_profile = []
    if r2:
        pos = collections.defaultdict(collections.Counter)
        m = 0
        for seq in fastq_seqs(r2):
            for i, b in enumerate(seq):
                pos[i][b] += 1
            m += 1
            if m >= r2_subsample:
                break
        for i in sorted(pos):
            c = pos[i]
            tot = sum(c.values())
            r2_profile.append({"pos": i, "n": tot,
                               "max_frac": max(c.values()) / tot,
                               "A": c["A"], "C": c["C"], "G": c["G"],
                               "T": c["T"], "N": c["N"]})

    qc = {"sample": name, "total_reads": n, "fate": dict(fate),
          "stagger": {STAGGERS.get(k, "other/unassigned"): v for k, v in stagger.items()},
          "read_len": dict(lens), "mapped_reads": sum(counts.values()),
          "distinct_guides_detected": len(counts),
          "unmatched_20mers_distinct": len(unmatched),
          "unmatched_top": unmatched.most_common(50),
          "r2_profile": r2_profile}
    with open(os.path.join(outdir, f"{name}_qc.json"), "w") as fh:
        json.dump(qc, fh)
    with open(os.path.join(outdir, f"{name}_counts.tsv"), "w") as fh:
        fh.write("guide_idx\tcount\n")
        for i, c in sorted(counts.items()):
            fh.write(f"{i}\t{c}\n")
    with open(os.path.join(outdir, f"{name}_unmatched.tsv"), "w") as fh:
        fh.write("seq\tcount\n")
        for s, c in unmatched.most_common():
            if c >= 5:
                fh.write(f"{s}\t{c}\n")
    print(f"[{name}] done: {n} reads, {sum(counts.values())} mapped", flush=True)
    return qc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", required=True)
    ap.add_argument("--r1", required=True)
    ap.add_argument("--r2", default=None)
    ap.add_argument("--key", required=True, help="CSV with Guide_Seq,Gene,Guide_Num")
    ap.add_argument("--outdir", default="results")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    guides = []
    with open(a.key) as fh:
        header = fh.readline().rstrip("\n").split(",")
        gi = header.index("Guide_Seq")
        for line in fh:
            guides.append(line.rstrip("\n").split(",")[gi].strip().upper())
    key_index = {g: i for i, g in enumerate(guides)}
    assert len(key_index) == len(guides), "duplicate guide sequences in key"
    count_sample(a.name, a.r1, a.r2, key_index, a.outdir)


if __name__ == "__main__":
    main()
