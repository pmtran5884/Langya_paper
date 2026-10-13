#!/usr/bin/env python3
"""Collapse the per-taxon outputs into the two Figure 5 / Supplementary Table 4 files.

Reads only what steps 01-04 already wrote under results/ - no re-extraction, no
re-alignment - and emits

    results/epcr_three_species.faa                the three recovered proteins
    results/supplementary_table_4.tsv             one row per taxon: accession,
                                                  locus, exon/CDS structure,
                                                  protein length, identity to
                                                  the miniprot query, QC flags
                                                  and the sequence itself

The three paper targets are the species with no annotated PROCR
(Blarina brevicauda, Crocidura gueldenstaedtii, Rattus flavipectus). Rattus
rattus is included as the control row: it has an annotation, and the pipeline's
recovered protein is byte-identical to RefSeq XP_032760446.1, which is what
makes it a check on the method rather than a result.
"""
import argparse
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# label in results/ -> (species, role, how the sequence was obtained)
TAXA = [
    ("blarina_brevicauda", "Blarina brevicauda", "target",
     "assembly, miniprot splice-aware protein-to-genome alignment"),
    ("crocidura_gueldenstaedtii", "Crocidura gueldenstaedtii", "target",
     "assembly, miniprot splice-aware protein-to-genome alignment"),
    ("SRR3948161", "Rattus flavipectus", "target",
     "WGS reads mapped to the R. rattus locus, consensus, miniprot"),
    ("rattus_rattus", "Rattus rattus", "control (annotated)",
     "assembly, miniprot splice-aware protein-to-genome alignment"),
]

COLS = ["species", "role", "label", "source", "accession", "contig",
        "cds_start", "cds_end", "strand", "cds_exons", "cds_bp", "protein_aa",
        "identity_to_query", "query_protein", "starts_with_M", "internal_stops",
        "frameshifts_reported", "n_terminal_trim", "mean_depth_over_gene",
        "protein_sequence"]


def read_kv(path):
    if not os.path.exists(path):
        return {}
    return dict(l.split("\t", 1) for l in open(path).read().splitlines() if "\t" in l)


def read_fasta_one(path):
    seq = []
    for line in open(path):
        if not line.startswith(">"):
            seq.append(line.strip())
    return "".join(seq)


def read_fasta(path):
    seqs, name = {}, None
    for line in open(path):
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = []
        elif name:
            seqs[name].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}


def cds_exons(gff):
    """Rank=1 CDS intervals from a miniprot GFF, 1-based inclusive."""
    ex = []
    if not os.path.exists(gff):
        return ex
    for line in open(gff):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) > 8 and f[2] == "CDS" and "Rank=1" in f[8]:
            ex.append((int(f[3]), int(f[4])))
    return sorted(ex)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--resdir", default=os.path.join(ROOT, "results"))
    a = ap.parse_args()
    res = a.resdir

    rows = []
    for label, species, role, source in TAXA:
        d = os.path.join(res, label)
        qc = read_kv(os.path.join(d, "qc.tsv"))
        prot = read_fasta_one(os.path.join(d, "protein.faa"))
        r = dict.fromkeys(COLS, "")
        r.update(species=species, role=role, label=label, source=source,
                 protein_sequence=prot, protein_aa=len(prot))
        if "assembly" in qc:                       # assembly route (step 01)
            r.update(accession=qc["assembly"], contig=qc["contig"],
                     cds_start=qc["cds_start"], cds_end=qc["cds_end"],
                     strand=qc["strand"], cds_exons=qc["cds_exons"],
                     cds_bp=qc["cds_bp"], identity_to_query=qc["identity_to_query"],
                     query_protein=qc["query"], starts_with_M=qc["starts_with_M"],
                     internal_stops=qc["internal_stops"],
                     frameshifts_reported=qc["frameshifts_reported"],
                     n_terminal_trim=qc["n_terminal_trim"])
        else:                                      # reads route (steps 02-03)
            ex = cds_exons(os.path.join(d, "miniprot.gff"))
            cds = read_fasta(os.path.join(res, "cds_all.fna")).get(label, "")
            r.update(accession=qc.get("accession", label),
                     mean_depth_over_gene=qc.get("gene_mean_depth", ""),
                     cds_bp=len(cds), cds_exons=len(ex),
                     query_protein="NP_001020904.1",
                     contig="reconstructed on the R. rattus NC_046158.1 locus",
                     cds_start=ex[0][0] if ex else "", cds_end=ex[-1][1] if ex else "",
                     strand="+ (consensus coordinates)",
                     internal_stops=prot.count("*"))
        rows.append(r)

    faa = os.path.join(res, "epcr_three_species.faa")
    with open(faa, "w") as fh:
        for r in rows:
            if r["role"] != "target":
                continue
            fh.write(f">{r['species'].replace(' ', '_')} EPCR|PROCR {r['accession']} "
                     f"{r['protein_aa']}aa\n{r['protein_sequence']}\n")

    tsv = os.path.join(res, "supplementary_table_4.tsv")
    with open(tsv, "w") as fh:
        fh.write("\t".join(COLS) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in COLS) + "\n")

    for r in rows:
        print(f"[05] {r['species']:<28} {r['protein_aa']} aa  {r['role']}")
    print(f"[05] wrote {faa}")
    print(f"[05] wrote {tsv}")


if __name__ == "__main__":
    main()
