# Figure 5 / Supplementary Table 4 — EPCR (PROCR) in three species with no annotation

Three target species have no annotated *PROCR*, so each EPCR sequence used in
the paper is a prediction made here, by one of two routes:

| species | route | source | result |
|---|---|---|---|
| *Blarina brevicauda* | assembly | GCA_051295465.2 | 242 aa, `JBKFBT020000153.1:87,596-91,813(+)`, 729 bp, 4 exons |
| *Crocidura gueldenstaedtii* | assembly | GCA_965212715.2 | 241 aa, `OZ243422.1:4,125,072-4,130,097(−)`, 726 bp, 4 exons |
| *Rattus flavipectus* | WGS reads | SRR3948161 (12.5× over the gene), SRR2917427 (32.2×) | 240 aa, 723 bp, 4 exons |
| *Rattus rattus* (control) | assembly | GCF_011064425.1 | 240 aa, identical to RefSeq XP_032760446.1 |

*R. rattus* is annotated and is included as the control on the method, not as a
result: the pipeline's protein for it is identical to the existing RefSeq
entry, and its locus also supplies the mapping scaffold for the reads route.

**Assembly route** (`01`): `miniprot` aligns a query EPCR protein to the raw
assembly with a splice model, so no GFF/GTF is needed; `extract_cds_from_gff.py`
lifts the rank-1 CDS, trims any N-terminal over-extension to the first in-frame
ATG, translates, and writes the QC row (frame, start codon, internal stops,
frameshifts, identity to query). Query is *Sorex araneus* XP_004612598.1 for the
shrews and *R. norvegicus* NP_001020904.1 for the rats.

**Reads route** (`02`, `03`): *R. flavipectus* has no assembly. Reads are
captured against the *R. rattus* PROCR locus ± 25 kb, competitively re-mapped to
the whole genome at MAPQ ≥ 20, called into a consensus (`--min-depth 3`,
`--min-BQ 15 --min-MQ 20`) in both a majority-base and an IUPAC flavour, and the
CDS is lifted from each by miniprot. This is reference-guided reconstruction,
not de novo assembly: it reports only positions where reads mapped, so a
*flavipectus*-specific insertion or a region too divergent to map would appear
as N rather than as a wrong base.

## Reproduce the tables and sequences from committed results (~20 s)

Everything under `results/` is committed, so the collection, alignment and
table steps re-run without downloading anything:

```bash
python scripts/04_align_and_report.py     # MAFFT alignments + identity tables + report.md
python scripts/05_supp_table4.py          # Supplementary Table 4 + the three-sequence FASTA
```

Needs only MAFFT on `PATH` plus Python 3 (standard library). Verified: the
regenerated `proteins_all.faa`, `cds_all.fna`, `protein_identity.tsv`,
`summary.tsv` and `codon_differences.tsv` are byte-identical to the committed
copies.

Deliverables:

* `results/epcr_three_species.faa` — the three predicted proteins, for Figure 5
* `results/supplementary_table_4.tsv` — one row per taxon: accession, locus,
  exon/CDS structure, protein length, identity to the miniprot query, QC flags,
  read depth where applicable, and the sequence
* `results/protein_identity.tsv` — all pairwise protein identities from the
  MAFFT alignment (the two shrews are 82.08 % to each other, 62.9-65.0 % to the
  rats; *R. flavipectus* and *R. rattus* are 100 %)
* `results/summary.tsv`, `codon_differences.tsv` — CDS-level comparison inside
  the *Rattus* group only, where the CDS are the same length (723 bp); the four
  *flavipectus* differences are all synonymous
* `results/report.md` — the same content as a readable report

## Re-run the extraction from scratch

```bash
./run_all.sh                       # everything, in order
scripts/01_extract_from_assembly.sh blarina_brevicauda        # one species
```

Cost and runtime, so nothing is a surprise:

* assembly route — one ~1-2 GB genome download per species, then minutes of
  miniprot. Both shrews plus the *R. rattus* control: well under an hour.
* reads route — 15 GB (SRR3948161) and 39 GB (SRR2917427) of `.sra` per run,
  then one streaming pass over all reads. ~1-2 h wall clock per run.

Large intermediates (genomes, `.sra`, BAMs) go to `work/`, which is gitignored;
small publishable outputs go to `results/`.

`02_fetch_sra.sh` does not use `prefetch`: sra-tools resolves accessions through
a name-resolver service that many networks block, after which even local files
fail with "Failed to call external services". It asks the NCBI SDL endpoint for
the delivery URL and md5, pulls the object with ranged parallel HTTP, verifies
the md5 against NCBI's manifest, and `03_reconstruct_locus.sh` then runs
sra-tools fully offline with its own `NCBI_SETTINGS`.

## Requirements

`miniprot`, `minimap2`, `samtools` (≥ 1.17, for `samtools consensus`),
`bcftools`, `sra-tools` (`fasterq-dump`), `mafft`, `curl`, Python 3. See
`../environment.yml`.

Species, modes, accessions and query proteins are in `config/targets.tsv`;
parameters (flank, threads, depth and quality floors, splice model) in
`config/params.sh`. Retargeting to another gene or taxon means editing only
those two files.

## Scope limits worth stating with these sequences

* *R. flavipectus* is reference-guided reconstruction (see above). Two
  independent runs/individuals at 12.5× and 32× agree, and all four CDS
  differences from *R. rattus* are synonymous, so the protein is identical.
* *Blarina brevicauda* rests on a single assembly with no independent
  replication and no annotation to check against — only internal QC (frame,
  start, stop, no internal stops, zero frameshifts reported). Nothing looks
  wrong — 86.9 % identity to the *Sorex araneus* query and 82.1 % to
  *Crocidura*, as the phylogeny predicts — but "no anomaly" is weaker than
  "checked".
* The assembly deposited under GCA_965212715 was relabelled between versions
  (`.1` = *Crocidura zimmermanni*, `.2` = *C. gueldenstaedtii*, internal name
  still `CroZim`), so cite the accession **with its version**, not the species
  name alone.
