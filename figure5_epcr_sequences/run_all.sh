#!/usr/bin/env bash
# EPCR recovery, end to end. Expect ~1-2 h wall clock per WGS run, dominated by
# the download (15-39 GB per run) and the single streaming pass over all reads;
# the assembly route is minutes per species after a ~1-2 GB genome download.
#
#   conda env create -f ../environment.yml && conda activate langya-paper
#   ./run_all.sh                 # everything
#   ./run_all.sh --tables-only   # steps 04-05 from the committed results/
set -euo pipefail
cd "$(dirname "$0")"
source config/params.sh

if [[ "${1:-}" == "--tables-only" ]]; then
  scripts/04_align_and_report.py
  scripts/05_supp_table4.py
  echo "[run_all] tables regenerated from the committed results/"
  exit 0
fi

# 1. the scaffold species: assembly + annotation, also supplies the bait
scripts/01_extract_from_assembly.sh "$SCAFFOLD_LABEL"

# 2. the two shrews: unannotated assemblies, miniprot finds the gene directly
for label in $GENOME_TARGETS; do
  scripts/01_extract_from_assembly.sh "$label"
done

# 3. the taxon with no assembly: reconstruct the locus from WGS reads
for acc in $RUNS; do
  scripts/02_fetch_sra.sh "$acc"
  scripts/03_reconstruct_locus.sh "$acc"
done

# 4. align everything, write the tables, then the paper deliverables
scripts/04_align_and_report.py
scripts/05_supp_table4.py
echo "[run_all] done -- see $RESULTS/report.md, $RESULTS/supplementary_table_4.tsv"
