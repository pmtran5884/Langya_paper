# Langya_paper

Code and inputs for the data-generating analyses in the Langya virus (LayV)
receptor paper. One directory per figure, each self-contained: scripts, the
input tables they need, and the outputs they produce, so every panel can be
re-run on its own.

| directory | panel | what it does | offline re-run |
|---|---|---|---|
| [`figure1A_parahenipavirus_map/`](figure1A_parahenipavirus_map) | Fig. 1A | world map of reported parahenipaviruses | yes, ~10 s (R) |
| [`figure1C_crispra_screen/`](figure1C_crispra_screen) | Fig. 1C | CRISPRa surfaceome screen for the LayV-G binding partner → volcano | yes, ~5 s (Python) |
| [`figure5_epcr_sequences/`](figure5_epcr_sequences) | Fig. 5, Supp. Table 4 | EPCR (PROCR) recovery in three species with no annotation | yes, ~20 s (Python + MAFFT) |

This is the minimal path to each panel. The fuller working repositories behind
two of them are `pmtran5884/Lay-G-CRISPRa` (screen QC, depth subsampling,
library-representation analysis, an independent reimplementation of the gene
scoring) and `pmtran5884/epcr-gene-recovery` (per-species audits, dN/dS, the
original standalone recovery scripts). Nothing here depends on either.

## Quick start — reproduce every output

```bash
./run_all.sh            # the three offline paths, in order
```

`run_all.sh` runs only what needs no large download and no network: the map,
the volcano from Apron's committed output, and the EPCR tables from the
committed per-species results. Each figure directory's README gives the full
from-scratch path as well (raw reads, genome downloads, a live portal
submission) with its runtime and disk cost.

## What is committed, and what is not

Committed, so a reviewer can check a number without re-running anything:

* the sgRNA library key and the per-guide count matrix for all three screen arms
* the exact `.chip` and log-fold-change files submitted to the Broad GPP portal,
  and the portal's unmodified gene table and volcano PDF
* every recovered EPCR protein and CDS with its per-species QC row, the MAFFT
  alignments, the identity tables and Supplementary Table 4
* the rendered Figure 1A panel

Not committed (`.gitignore`): raw FASTQ, `.sra` files, genome assemblies, BAMs
and other large intermediates. Raw sequencing data goes to SRA at publication.
`figure1C_crispra_screen/example_data/` holds a 30,000-read subsample of two
screen arms purely so the counting step is runnable; it is too shallow to
reproduce the paper's gene scores.

## Verification already done

These are checks on this repository, not on the original analyses:

* `figure1C .../02_make_apron_input.py` reproduces the submitted
  `layg_crispra.chip` and `layg_crispra_input.txt` **byte for byte** from the
  committed count matrix.
* `figure1C .../04_volcano_apron.py` reproduces the paper's gene-score table
  exactly (max absolute difference 0 against the table in `Lay-G-CRISPRa`), with
  PROCR at Apron score 9.145947, p = 7.15 × 10⁻¹⁰, BH-FDR = 2.22 × 10⁻⁶ — the
  only gene below FDR 0.1.
* `figure1C .../01b_build_count_matrix.py` round-trips the committed count
  matrix exactly.
* `figure5 .../04_align_and_report.py` regenerates `proteins_all.faa`,
  `cds_all.fna`, `protein_identity.tsv`, `summary.tsv` and
  `codon_differences.tsv` byte-identically to the committed copies.
* `figure1A .../make_world_map.R` runs clean and reproduces the committed panel.

## Environment

```bash
conda env create -f environment.yml
conda activate langya-paper
```

Figure 1A needs R with `ggplot2`, `sf`, `dplyr`, `rnaturalearth`,
`rnaturalearthdata` (not in the conda file; see that directory's README).
Figure 1C needs Python with `numpy`, `pandas`, `matplotlib`, `requests`.
Figure 5 needs `mafft` for the offline path, plus `miniprot`, `minimap2`,
`samtools`, `bcftools` and `sra-tools` for a full re-run.

## Key results reproduced here

* **Fig. 1C** — PROCR (EPCR) is the only gene below FDR 0.1 in the CRISPRa
  screen (p = 7.15 × 10⁻¹⁰), all four guides enriched, read share rising from
  0.005 % in the unselected library to 2.04 % in the sorted arm.
* **Fig. 5 / Supp. Table 4** — EPCR predicted for three species with no
  annotated *PROCR*: *Blarina brevicauda* (242 aa), *Crocidura gueldenstaedtii*
  (241 aa) and *Rattus flavipectus* (240 aa), with *Rattus rattus* as the
  control whose predicted protein matches RefSeq XP_032760446.1.
