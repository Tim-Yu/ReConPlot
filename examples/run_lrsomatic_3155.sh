#!/usr/bin/env bash
# Worked example: lrsomatic output for P215003155 (Severus SVs + ASCAT CN).
#
#   bash examples/run_lrsomatic_3155.sh
#
# Produces, under $OUT:
#   1. per_chromosome/       one figure per chromosome, full somatic callset
#   2. genome_wide/          overview strip
#   3. focus_cluster/        chr16+chr17 (Severus cluster severus_0) with a BAF
#                            track and ASCAT's unrounded raw segmentation
#   4. rearrangements_only/  the same view with sub-kb and VNTR calls removed,
#                            for like-for-like comparison with SAVANA
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
IN=/mnt/dish/0.tempshare_ToDel/ONT3155bam/lrsomatic/output/P215003155
OUT=/mnt/scratch/BYU/6.savana/data/ReConPlot_output/P215003155_lrsomatic

"$HERE/run_reconplot.R" --source lrsomatic --input "$IN" --outdir "$OUT/per_chromosome" \
  --regions all --layout separate --format pdf,png --write-tables

"$HERE/run_reconplot.R" --source lrsomatic --input "$IN" --outdir "$OUT/genome_wide" \
  --regions all --layout together --width 34 --height 4 --format pdf,png \
  --extra "size_text=4,size_chr_labels=5,scale_ticks=50000000"

# chr16+chr17 carry Severus cluster severus_0 (191 + 48 breakpoints) and the
# 36 chr16<->chr17 translocations, so they belong on one panel.
"$HERE/run_reconplot.R" --source lrsomatic --input "$IN" --outdir "$OUT/focus_cluster" \
  --regions "chr16,chr17" --layout together --baf-track \
  --ascat-cn-mode raw-fractional --max-cn 6 \
  --genes TP53,ERBB2,BRCA1,NF1 --format pdf,png

# Severus runs at --min-support 3 and reports small VNTR indels next to
# rearrangements; drop those to compare against the SAVANA figures.
"$HERE/run_reconplot.R" --source lrsomatic --input "$IN" --outdir "$OUT/rearrangements_only" \
  --regions "chr16,chr17" --layout together \
  --min-svlen 1000 --exclude-vntr \
  --ascat-cn-mode raw-fractional --max-cn 6 --format pdf,png

echo "Wrote figures to $OUT"
