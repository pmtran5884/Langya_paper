#!/usr/bin/env bash
# Reproduce every committed output, using only what is in this repository:
# no large downloads, no network, no raw reads. Expect well under a minute,
# plus whatever R takes to load the spatial stack.
#
#   conda env create -f environment.yml && conda activate langya-paper
#   ./run_all.sh
#
# Each figure directory's README gives the full from-scratch path as well.
set -euo pipefail
cd "$(dirname "$0")"
ROOT="$(pwd)"

echo "=== Figure 1A: parahenipavirus world map"
if command -v Rscript >/dev/null 2>&1; then
  mkdir -p figure1A_parahenipavirus_map/output
  ( cd figure1A_parahenipavirus_map/output && Rscript ../make_world_map.R )
  echo "--> figure1A_parahenipavirus_map/output/Parahenipavirus_world_map_saturated.{png,pdf}"
else
  echo "    SKIPPED: Rscript not on PATH (see figure1A_parahenipavirus_map/README.md)"
fi

echo
echo "=== Figure 1C: CRISPRa screen volcano, from Apron's committed output"
( cd figure1C_crispra_screen && python scripts/04_volcano_apron.py )

echo
echo "=== Figure 5 / Supp. Table 4: EPCR tables, from the committed per-species results"
if command -v mafft >/dev/null 2>&1; then
  ( cd figure5_epcr_sequences && ./run_all.sh --tables-only )
else
  echo "    alignment step SKIPPED: mafft not on PATH; writing the table only"
  ( cd figure5_epcr_sequences && python scripts/05_supp_table4.py )
fi

echo
echo "=== done. Outputs:"
echo "  figure1A_parahenipavirus_map/output/Parahenipavirus_world_map_saturated.{png,pdf}"
echo "  figure1C_crispra_screen/figures/volcano_apron_LayG1st.{png,pdf,svg}"
echo "  figure1C_crispra_screen/results/apron_gene_scores_LayG1st_vs_Mock.csv"
echo "  figure5_epcr_sequences/results/epcr_three_species.faa"
echo "  figure5_epcr_sequences/results/supplementary_table_4.tsv"
echo "  figure5_epcr_sequences/results/report.md"
