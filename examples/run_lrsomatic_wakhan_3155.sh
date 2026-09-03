#!/usr/bin/env bash
# Worked example: lrsomatic output for P215003155 (Severus SVs + Wakhan CN).
#
#   bash examples/run_lrsomatic_wakhan_3155.sh
#
# Produces, under $OUT:
#   1. per_chromosome/       one figure per chromosome, full somatic callset
#   2. genome_wide/          overview strip using Wakhan allele-specific CN
#   3. focus_cluster/        chr16+chr17 (Severus cluster severus_0)
#   4. rearrangements_only/  the same view with sub-kb and VNTR calls removed
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
IN=/mnt/dish/0.tempshare_ToDel/ONT3155bam/lrsomatic/output/P215003155
OUT=/mnt/scratch/BYU/6.savana/data/ReConPlot_output/P215003155_lrsomatic_wakhan

"$HERE/run_reconplot.R" --cn-source wakhan --sv-source severus --input "$IN" --outdir "$OUT/per_chromosome" \
  --regions all --layout separate --format pdf,png --write-tables

"$HERE/run_reconplot.R" --cn-source wakhan --sv-source severus --input "$IN" --outdir "$OUT/genome_wide" \
  --regions all --layout together --width 34 --height 4 --format pdf,png \
  --extra "size_text=4,size_chr_labels=5,scale_ticks=50000000"

"$HERE/run_reconplot.R" --cn-source wakhan --sv-source severus --input "$IN" --outdir "$OUT/focus_cluster" \
  --regions "chr16,chr17" --layout together --max-cn 8 \
  --genes TP53,ERBB2,BRCA1,NF1 --format pdf,png

"$HERE/run_reconplot.R" --cn-source wakhan --sv-source severus --input "$IN" --outdir "$OUT/rearrangements_only" \
  --regions "chr16,chr17" --layout together \
  --min-svlen 1000 --exclude-vntr --max-cn 8 --format pdf,png

echo "Wrote figures to $OUT"
