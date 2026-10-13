# Figure 1C — CRISPRa surfaceome screen for the LayV-G binding partner

A375 cells carrying the CRISPRa surfaceome library were stained with LayV-G and
sorted; the sorted arm (`LayG1st`) is compared against the unselected library
(`Mock`). Gene scoring is done by the Broad GPP pooled-screen tool (Apron), not
by us — the scripts here prepare its input and render its output.

Four steps, each runnable on its own:

| step | script | in | out |
|---|---|---|---|
| 01 | `scripts/01_count_sgrna.py` | R1 FASTQ + library key | per-sample `*_counts.tsv`, `*_qc.json` |
| 01b | `scripts/01b_build_count_matrix.py` | per-sample counts | `sgrna_count_matrix.csv` |
| 02 | `scripts/02_make_apron_input.py` | count matrix | `.chip` + log-fold-change `.txt` for the portal |
| 03 | `scripts/03_submit_apron.py` | those two files | Apron's gene table |
| 04 | `scripts/04_volcano_apron.py` | Apron's gene table | gene-score CSV + **Figure 1C** |

## Reproduce the figure (offline, ~5 s)

`apron_submission/` holds the exact files submitted for the paper and the
portal's unmodified output, so step 04 needs no network and no raw reads:

```bash
python scripts/04_volcano_apron.py
```

→ `figures/volcano_apron_LayG1st.{png,pdf,svg}` (210 × 130 mm, Arial 10 pt) and
`results/apron_gene_scores_LayG1st_vs_Mock.csv` (all 6,138 scored genes).

Expected console line:

```
LayG1st: 6,138 genes (249 control) | enriched FDR<0.1: 1 | PROCR score 9.145947 p=7.15e-10 FDR=2.22e-06
```

PROCR (EPCR) is the only gene below FDR 0.1 in either direction.

## Reproduce the analysis from the count matrix (~10 s, needs network)

```bash
python scripts/02_make_apron_input.py --outdir my_submission
python scripts/03_submit_apron.py --chip my_submission/layg_crispra.chip \
       --data my_submission/layg_crispra_input.txt --outdir my_apron_out
python scripts/04_volcano_apron.py --apron-out my_apron_out
```

Step 02 is checked: run as above it reproduces `apron_submission/layg_crispra.chip`
and `layg_crispra_input.txt` **byte for byte**. Step 03 drives the public web
form; the same two files can instead be uploaded by hand at
<https://portals.broadinstitute.org/gpp/public/analysis-tools/crispr-gene-scoring>
with the paper's settings (printed by step 02, and listed below).

## Reproduce the counts from reads

Raw FASTQs are not in git (they go to SRA at publication). `example_data/` holds
a 30,000-read subsample of two arms so the counting step is runnable:

```bash
python scripts/01_count_sgrna.py --name Mock    --r1 example_data/example_mock_R1.fastq.gz \
       --key data/formatted_sgrna_key.csv --outdir counts
python scripts/01_count_sgrna.py --name LayG1st --r1 example_data/example_layg1st_R1.fastq.gz \
       --key data/formatted_sgrna_key.csv --outdir counts
python scripts/01b_build_count_matrix.py --countdir counts --samples Mock LayG1st \
       --out counts/example_count_matrix.csv
```

A subsample this small will not reproduce the paper's gene scores — the full
arms are 12.4 M (Mock), 4.4 M (LayG1st) and 3.6 M (LayG2nd) mapped reads. The
paper's counts are committed as `data/sgrna_count_matrix.csv`, which is where
steps 02-04 start.

R1 architecture the counter anchors on (the sgRNA start is not at a fixed
offset, because the four px52 forward primers carry a 0-3 nt stagger):

```
[adapter] [0-3 nt stagger] CTTGTGGAAAGGACGAAACACCG [20 nt sgRNA] GTTTAAGAGCTATGCTG
```

R2 is uninformative (>95 % invariant at every position) and is used only for the
QC base-composition profile.

## Analysis choices that change the ranking

* **Normalisation is ours, not Apron's.** The portal consumes log-fold changes
  and never sees counts: `lognorm = log2(CPM + 1)`, `LFC = lognorm(arm) −
  lognorm(Mock)`, totals over all 24,367 guides.
* **1,039 guides with zero Mock reads are omitted from the data file** (they stay
  in the `.chip`). Their fold change is undefined, and Apron has no
  minimum-reference-reads option; leaving them in ranks library-absent genes
  above genuinely depleted ones.
* **Non-targeting guides are renamed `NONTARGET_n` → `NO_SITE_n`**, the portal's
  convention. Without it Apron scores each as a one-guide gene instead of
  pooling the 1,000 into 249 four-guide control pseudo-genes.
* **Portal settings:** Hypergeometric, directionality Both, top 100 % of guides,
  dummy gene set size 4 (library mean is 3.97 guides/gene), guides per gene
  min 2 / max 4, no-site controls displayed.
* **The BH-FDR is an addition.** The GPP hypergeometric output carries no
  corrected p-value; step 04 computes BH within each direction and names the
  column `fdr_BH_added` to keep that visible.

## Files

* `data/formatted_sgrna_key.csv` — library key, 24,367 guides (`Guide_Seq`, `Gene`, `Guide_Num`)
* `data/sgrna_count_matrix.csv` — per-guide counts for Mock, LayG1st, LayG2nd
* `apron_submission/layg_crispra.chip`, `layg_crispra_input.txt` — exactly what was submitted
* `apron_submission/LayG1st_vs_Mock_100_26-10-03-09-05.txt` — Apron's unmodified gene table
* `apron_submission/LayG1st_vs_Mock.pdf` — Apron's own volcano, for comparison
* `scripts/figure_style.py` — figure defaults shared by step 04
