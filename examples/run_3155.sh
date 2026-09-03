#!/usr/bin/env bash
# Worked example: SAVANA output for P215003155_tumor.
#
#   bash examples/run_3155.sh
#
# Produces, under $OUT:
#   1. one figure per chromosome (PDF + PNG)
#   2. a genome-wide overview strip
#   3. a focused chr17 chromothripsis panel with a BAF track and gene labels
#   4. the harmonised CN/SV tables used for all of the above
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
IN=/mnt/scratch/BYU/6.savana/data/3155
OUT=/mnt/scratch/BYU/6.savana/data/ReConPlot_output/P215003155_tumor

"$HERE/run_reconplot.R" --cn-source savana --sv-source savana --input "$IN" --outdir "$OUT/per_chromosome" \
  --regions all --layout separate --format pdf,png --write-tables

"$HERE/run_reconplot.R" --cn-source savana --sv-source savana --input "$IN" --outdir "$OUT/genome_wide" \
  --regions all --layout together --width 34 --height 4 --format pdf,png \
  --extra "size_text=4,size_chr_labels=5,scale_ticks=50000000"

"$HERE/run_reconplot.R" --cn-source savana --sv-source savana --input "$IN" --outdir "$OUT/focus" \
  --regions "chr17" --layout together --baf-track \
  --genes TP53,ERBB2,BRCA1,NF1 --format pdf,png

echo "Wrote figures to $OUT"
