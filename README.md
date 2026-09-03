# ReConPlot driver scripts

Wrapper around [ReConPlot](https://github.com/cortes-ciriano-lab/ReConPlot) that
turns SV/CN caller output into publication-style rearrangement plots.

ReConPlot itself is a single function; almost all of the work is getting the two
input tables into the exact shape it expects. That reshaping is what lives here,
split so that supporting a new caller means writing one small parser file.

```
<caller output>  --[ R/parsers/<caller>.R ]-->  {cn, sv}
                 --[ R/utils_common.R validate ]-->  ReConPlot::ReConPlot()
```

A parser may supply both tables (SAVANA) or just one (ASCAT supplies CN,
Severus supplies SVs). Components compose on the command line, so any CN caller
pairs with any SV caller:

```bash
--source lrsomatic                          # convenience wrapper
--cn-source ascat --sv-source severus       # the same thing, explicit
--cn-source savana --sv-source severus      # mix and match
```

## Layout

| Path | Role |
|---|---|
| `run_reconplot.R` | CLI entry point; caller-agnostic |
| `R/utils_common.R` | Data contracts, chromosome handling, region parsing |
| `R/plotting.R` | Sizing, titles, file naming, per-panel error handling |
| `R/utils_vcf.R` | Minimal VCF reading shared by the SV parsers |
| `R/parsers/registry.R` | `register_parser()` plug-in table + source composition |
| `R/parsers/savana.R` | SAVANA (CN + SVs) |
| `R/parsers/ascat.R` | ASCAT allele-specific copy number (CN only) |
| `R/parsers/severus.R` | Severus long-read somatic SVs (SVs only) |
| `R/parsers/lrsomatic.R` | lrsomatic pipeline = ASCAT + Severus |
| `R/parsers/generic.R` | Any caller, via explicit files + column auto-mapping |
| `examples/run_3155.sh` | Worked example: SAVANA on P215003155 |
| `examples/run_lrsomatic_3155.sh` | Worked example: lrsomatic on P215003155 |

## The two data contracts

Everything a parser returns is validated against these before plotting:

* **CN** — `chr`, `start`, `end`, `copyNumber`, `minorAlleleCopyNumber`
* **SV** — `chr1`, `pos1`, `chr2`, `pos2`, `strands`, where `strands` is one of
  `+-` (DEL-like), `-+` (DUP-like), `++` / `--` (INV-like), `TRA`, `INS`, `SBE`

Validation normalises chromosome names, drops alt/random/decoy contigs
(ReConPlot rejects them outright), coerces types, and by default recolours
inter-chromosomal junctions as `TRA`.

## Usage

```bash
./run_reconplot.R --list-sources

# whole-genome, one PDF+PNG per chromosome
./run_reconplot.R --source savana \
  --input  /mnt/scratch/BYU/6.savana/data/3155 \
  --outdir /mnt/scratch/BYU/6.savana/data/ReConPlot_output/P215003155_tumor

# a focused multi-panel figure with gene labels
./run_reconplot.R --source savana --input DIR --outdir OUT \
  --regions "chr8,chr17:30000000-50000000" --layout together \
  --genes MYC,TP53,ERBB2
```

### Key options

| Option | Meaning |
|---|---|
| `--source` | Parser to use (`savana`, `lrsomatic`, `generic`) |
| `--cn-source` / `--sv-source` | Compose a CN parser with an SV parser |
| `--input` | Caller output directory; files are discovered inside |
| `--cn-input` / `--sv-input` | Per-component directories when they differ |
| `--sample` | Sample prefix, if the directory holds more than one |
| `--cn-file` / `--sv-file` | Bypass discovery with explicit paths |
| `--sv-format` | SAVANA SVs from `bedpe` (default) or `vcf` |
| `--min-support` | Drop SVs below this tumour read support |
| `--min-svlen` | Drop short intra-chromosomal SVs (translocations are kept) |
| `--exclude-vntr` | Drop SVs flagged inside a VNTR (Severus) |
| `--cluster-id` | Keep only SVs in the named Severus cluster(s) |
| `--clustered-only` | Keep only SVs Severus assigned to some cluster |
| `--ascat-cn-mode` | `segments` (default), `raw`, or `raw-fractional` |
| `--regions` | `all` (default), `chr8`, `chr8:120000000-130000000`, comma-separated |
| `--regions-file` | BED file of regions instead of `--regions` |
| `--layout` | `separate` (one plot per region), `together`, `both` |
| `--genes` | Comma-separated HUGO symbols to label |
| `--baf-track` | Add a het-SNP BAF panel below the plot (SAVANA) |
| `--max-cn` | Copy number axis ceiling (default 8) |
| `--format` | `pdf`, `png`, or `pdf,png` (default) |
| `--write-tables` | Also dump the harmonised CN/SV tables as TSV |
| `--extra` | Escape hatch: `key=value,...` forwarded to `ReConPlot()` |

Note that `--regions all --layout separate` writes 24 figures; `--layout
together` on many chromosomes produces one very wide panel strip, so pass
`--width` explicitly there.

## Adding another caller

Create `R/parsers/<name>.R` with a function taking the parsed CLI options and
returning the two tables, then register it:

```r
parse_mycaller <- function(args) {
  list(cn   = <data.frame chr,start,end,copyNumber,minorAlleleCopyNumber>,
       sv   = <data.frame chr1,pos1,chr2,pos2,strands>,
       meta = list(sample = args$sample),      # optional; used in the title
       annotation_fn = NULL)                   # optional function(regions) -> chr,pos,y
}
register_parser("mycaller", parse_mycaller, "one-line description")
```

For a caller that only does one half, declare it and return only that table:

```r
register_parser("mycncaller", parse_mycncaller, "...", provides = "cn")
register_parser("mysvcaller", parse_mysvcaller, "...", provides = "sv")
```

`load_parsers()` picks the file up automatically — no other file changes.
Return raw-ish columns; type coercion, contig filtering and orientation
normalisation happen centrally in `validate_cn()` / `validate_sv()`.

Reusable helpers already available: `normalize_chrom()`, `bedpe_point()`
(handles both `start == end == POS` and spec BEDPE `start == POS-1`), the
synonym mapper `.syn()` in `generic.R`, and the VCF helpers in `R/utils_vcf.R`
— `vcf_read_records()`, `vcf_info_get()`, `vcf_info_flag()`, `vcf_alt_mate()`
(decodes `t[chr:pos[` breakend ALTs), `vcf_format_get()` and
`dedupe_breakend_pairs()` (collapses mate-paired BND records to one junction).

## SAVANA specifics

Files discovered under `--input`:

| Pattern | Use |
|---|---|
| `*_segmented_absolute_copy_number.tsv` | CN track (requires SAVANA's CN fit step) |
| `*.classified.somatic.bedpe` | SV track (default) |
| `*.classified.somatic.vcf` | SV track with `--sv-format vcf` |
| `*_fitted_purity_ploidy.tsv` | purity/ploidy shown in the title |
| `*_allele_counts_hetSNPs.bed` | BAF panel with `--baf-track` |

SAVANA's `BP_NOTATION` already uses ReConPlot's orientation convention, so
orientations pass through unchanged; only `<INS>`/`<SBE>` are renamed. In the
BEDPE, orientation and read support are packed into the name column
(`ID_39590|1608bp|TUMOUR_12|--`) and are unpacked by `savana_split_bedpe_name()`.
The BEDPE and VCF paths were checked to yield identical junction sets.

## lrsomatic specifics (Severus + ASCAT)

`--source lrsomatic --input <outdir>/<sample>` resolves both halves itself.
`--input` also accepts the pipeline `output/` directory or the run directory
above it when there is only one sample.

### Which files to use, and why

| Tool | File used | Why |
|---|---|---|
| Severus | `variants/severus/somatic_SVs/severus_somatic.vcf.gz` | the complete somatic callset, with orientation in `INFO/STRANDS` |
| ASCAT | `ascat/<sample>.segments.txt` | ASCAT's fitted integer allele-specific calls (`ascat.output$segments`) |
| ASCAT | `ascat/<sample>.purityploidy.txt` | purity (aberrant cell fraction) and ploidy for the title |
| ASCAT | `ascat/<sample>.tumour_tumourBAF.txt` | het-SNP BAF panel with `--baf-track` |

Files deliberately **not** used as the primary input:

* `somatic_SVs/breakpoint_clusters.tsv` — lists only the junctions Severus
  assigned to a complex-rearrangement cluster. On P215003155 that is 213 of the
  330 junctions in the VCF (every cluster junction also appears in the VCF; 117
  VCF junctions belong to no cluster). It is a filter/annotation layer, not the
  callset. Cluster membership is in the VCF as `INFO/CLUSTERID` anyway, so
  `--cluster-id severus_0` filters without reading the TSV. Its companion
  `breakpoint_clusters_list.tsv` is read only to log the largest clusters.
* `all_SVs/severus_all.vcf.gz` — germline included; available via
  `--severus-all` if you want it.
* `ascat/<sample>.cnvs.txt` — byte-for-byte `segments.txt` minus the sample
  column (`ascat.output$segments[2:6]` in the module), kept only as a fallback.
* `ascat/<sample>.segments_raw.txt` — the pre-fit ASPCF segmentation. Reachable
  via `--ascat-cn-mode raw` (rounded) or `raw-fractional` (unrounded
  `nAraw`/`nBraw`); see the caveat below.

### Orientation conventions

Severus `INFO/STRANDS` uses the same convention as ReConPlot and SAVANA
(DEL `+-`, DUP `-+`, INV `++`/`--`), so orientations pass straight through.
This was checked on junctions both callers report, e.g. `chr4:89170-90778`,
which Severus writes as `STRANDS=--` and SAVANA as `--`. The sign prefixes in
`breakpoint_clusters.tsv` (`-chr17:47563030|-chr17:47608170`) concatenate to
the same value. ASCAT writes chromosomes without a `chr` prefix; that is
normalised on the way in.

### Two caveats worth knowing

* **Severus reports small VNTR indels.** lrsomatic runs it at
  `--min-support 3`, so the somatic VCF mixes rearrangements with sub-kb
  indels: 99 junctions under 1 kb, versus 19 from SAVANA on the same sample.
  Nothing is filtered by default. Use `--min-svlen 1000 --exclude-vntr` for a
  rearrangement-only view.
* **ASCAT segmentation is much coarser than SAVANA's.** 90 segments genome-wide
  against SAVANA's 553; on chr17 it is 3 against 54, so the fitted calls flatten
  the chromothripsis oscillation to a 2 -> 4 step. `--ascat-cn-mode
  raw-fractional` recovers some structure (9 segments on chr17 with unrounded
  values) but cannot match a read-depth caller. This is a property of ASCAT on
  G1000 SNP loci, not of these scripts.

## Requirements

R 4.4 with `ReConPlot`, `ggplot2`, `data.table`, `optparse`, `cowplot`,
`GenomicRanges`, `ggforce` — all present on this machine. Gzipped VCFs are read
through R's own `gzfile()`, so no tabix/bcftools is needed. Gene labels use
coordinates bundled with ReConPlot, so no network access is needed.

## Known ReConPlot quirks handled here

* Non-primary contigs in either table abort the call — filtered out up front.
* A `strands` value outside the known vocabulary aborts the call — filtered out.
* Repeating a chromosome in one multi-panel selection fails with
  "factor level [n] is duplicated" — repeated chromosomes are merged into their
  spanning window.
* Harmless warnings come from ggplot2 deprecations inside the package and from
  `min()`/`max()` over insertions on chromosomes outside the current selection.
