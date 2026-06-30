# Maize GWAS and LD Decay Plotting Scripts

This repository contains R scripts used to generate publication-ready visualizations and summary tables from TASSEL GWAS and linkage disequilibrium (LD) output files for maize association mapping analyses.

The scripts were developed for maize GWAS datasets exported from TASSEL and produce:

* Manhattan plots
* QQ plots with confidence envelopes
* Cleaned marker summary tables
* Top GWAS hit tables
* LD decay curves
* LD decay distance tables
* LD decay plots with or without raw pairwise LD point clouds

## Repository contents

```text
.
├── GWAS_graph.R
├── GWAS_LD_graph.R
├── GWAS_LD_scatter_graph.R
└── README.md
```

## Script overview

### `GWAS_graph.R`

Generates Manhattan and QQ plots from TASSEL GWAS output files.

Main outputs:

```text
plots/
├── <dataset>_manhattan.png
├── <dataset>_manhattan.svg
├── <dataset>_qq.png
└── <dataset>_qq.svg

tables/
├── <dataset>_cleaned_chromosome_markers.csv
├── <dataset>_top_hits.csv
├── <dataset>_chromosome_marker_summary.csv
└── CTI_summary_statistics.csv
```

The script keeps only chromosome-assigned, valid markers and removes markers with missing chromosome, position, marker name, or invalid p-values.

The Manhattan plots include:

* chromosome-wise coloring,
* cumulative genomic position on the x-axis,
* `-log10(p)` on the y-axis,
* Bonferroni threshold line,
* suggestive threshold line.

The QQ plots include:

* observed vs expected `-log10(p)` values,
* diagonal null expectation line,
* genomic inflation factor (`lambda GC`),
* null-distribution confidence envelope based on the beta distribution of ordered p-values.

### `GWAS_LD_graph.R`

Generates LD decay plots and LD decay tables from TASSEL LD output files.

This version plots chromosome-wise fitted or smoothed LD decay curves without displaying the full raw LD point cloud.

Implemented LD decay methods:

* `rolling_mean`: binned mean r² followed by moving-average smoothing,
* `loess`: LOESS curve fitted to binned mean r² values,
* `raw_bin_mean`: unsmoothed binned mean r²,
* `hill_weir`: Hill-Weir / Remington expected r² model.

Main outputs:

```text
LD_decay_binned_values.csv
LD_decay_<method>.png
LD_decay_<method>.svg
LD_decay_<method>_curve_values.csv
LD_decay_<method>_table.csv
LD_decay_<method>_decay_diagnostics.csv
```

### `GWAS_LD_scatter_graph.R`

Generates LD decay plots and LD decay tables from TASSEL LD output files, similarly to `GWAS_LD_graph.R`, but additionally displays a grey raw pairwise LD point cloud behind the fitted chromosome-wise LD decay curves.

This is useful for producing Remington-style LD decay figures where each grey point represents one raw pairwise TASSEL r² value and the colored curves show chromosome-wise fitted or smoothed LD decay trends.

The raw point cloud is used only for visualization. Decay distances and fitted curves are calculated from the selected LD decay method.

## Requirements

The scripts are written in R and require the following R packages.

For `GWAS_graph.R`:

```r
tidyverse
readr
ggplot2
svglite
ragg
scales
```

For `GWAS_LD_graph.R` and `GWAS_LD_scatter_graph.R`:

```r
data.table
dplyr
ggplot2
zoo
readr
tidyr
```

The `GWAS_graph.R` script automatically checks for missing packages and installs them from CRAN. The LD scripts assume that the required packages are already installed.

To install all dependencies manually:

```r
install.packages(c(
  "tidyverse",
  "readr",
  "ggplot2",
  "svglite",
  "ragg",
  "scales",
  "data.table",
  "dplyr",
  "zoo",
  "tidyr"
))
```

## Input data

### GWAS input format

`GWAS_graph.R` expects TASSEL GWAS output files in CSV-like format.

Required columns:

```text
Trait
Marker
Chr
Pos
p
```

Optional columns that are retained in the top-hit table if present:

```text
F
MarkerR2
df
errordf
Genetic Var
Residual Var
-2LnLikelihood
```

The script applies the following filters:

* `Marker` must not be missing or `None`,
* `Chr` must not be missing,
* `Pos` must be numeric and non-negative,
* `p` must be numeric and within `0 < p <= 1`,
* only chromosomes listed in `chr_order_letters` are retained,
* duplicate marker records are removed based on `Chr`, `Pos`, and `Marker`.

By default, the script expects chromosome labels `A` to `J`, corresponding to 10 maize chromosomes. The axis labels can optionally be changed to `1` to `10`.

### LD input format

`GWAS_LD_graph.R` and `GWAS_LD_scatter_graph.R` expect TASSEL LD export files.

Required columns:

```text
Locus1
Locus2
Dist_bp
R^2
```

Optional columns:

```text
pDiseq
N
```

The LD scripts keep only marker pairs where:

* `Locus1` and `Locus2` are on the same chromosome,
* `Dist_bp` is valid,
* `R^2` is valid and within `0 <= R^2 <= 1`,
* distance is within the user-defined maximum distance,
* optional `pDiseq` and `N` filters are satisfied if enabled.

The `hill_weir` method requires a valid `N` column.

## How to use

### 1. Edit input file paths

Each script contains a `USER SETTINGS` section near the top. Before running, edit the input paths to match your local files.

For example, in `GWAS_graph.R`:

```r
input_files <- c(
  SPAD2018 = "/path/to/SPAD_2018.csv.txt",
  SPAD2022 = "/path/to/SPAD_2022.csv.txt"
)
```

The names on the left side, such as `SPAD2018` and `SPAD2022`, are used as dataset labels in output file names and plot titles.

For the LD scripts, edit:

```r
input_csv <- "/path/to/LD_maize_TASSEL.csv.txt"
out_prefix <- "LD_decay"
```

### 2. Run the GWAS plotting script

From the command line:

```bash
Rscript GWAS_graph.R
```

This creates:

```text
plots/
tables/
```

and writes Manhattan plots, QQ plots, cleaned marker tables, top-hit tables, chromosome marker summaries, and an overall summary statistics table.

### 3. Run the LD decay script without raw point cloud

```bash
Rscript GWAS_LD_graph.R
```

This generates LD decay curves and tables using the method selected in the script.

### 4. Run the LD decay script with raw point cloud

```bash
Rscript GWAS_LD_scatter_graph.R
```

This generates LD decay curves with a grey background point cloud showing raw pairwise LD values.

## Important user settings

### GWAS script settings

The most important settings in `GWAS_graph.R` are:

```r
input_files
output_plot_dir
output_table_dir
rename_chr_to_numbers
qq_confidence_level
manhattan_width
manhattan_height
qq_width
qq_height
png_dpi
n_top_hits
```

#### Chromosome labels

By default:

```r
rename_chr_to_numbers <- FALSE
```

This keeps TASSEL chromosome labels `A` to `J`.

To display chromosomes as `1` to `10` on the Manhattan plot x-axis:

```r
rename_chr_to_numbers <- TRUE
```

#### QQ confidence envelope

The QQ plot confidence envelope is controlled by:

```r
qq_confidence_level <- 0.95
```

This produces a 95% confidence envelope under the null distribution of ordered p-values.

#### Top GWAS hits

The number of top markers exported to the top-hit table is controlled by:

```r
n_top_hits <- 50
```

### LD decay script settings

The most important LD settings are:

```r
max_dist_kb <- 2000
bin_size_kb <- 1
rolling_window <- 7
loess_span <- 0.25
loess_degree <- 2
thresholds <- c(0.2, 0.1)
decay_method <- "hill_weir"
compare_all_methods <- FALSE
```

#### Maximum distance

```r
max_dist_kb <- 2000
```

Only marker pairs up to this distance are included in plotting and LD decay estimation.

#### LD decay thresholds

```r
thresholds <- c(0.2, 0.1)
```

The scripts estimate the physical distance where the fitted or smoothed LD curve crosses each selected r² threshold.

For example, if `thresholds <- c(0.2, 0.1)`, the output table will contain columns similar to:

```text
LD_decay_kb_r2_0.2
LD_decay_kb_r2_0.1
```

#### LD decay method

Choose one method:

```r
decay_method <- "hill_weir"
```

Available options:

```r
"rolling_mean"
"loess"
"raw_bin_mean"
"hill_weir"
```

To run all methods:

```r
compare_all_methods <- TRUE
```

#### Optional LD filtering

To keep only significant LD pairs based on TASSEL `pDiseq`:

```r
use_pDiseq_filter <- TRUE
pDiseq_max <- 0.001
```

To require a minimum sample size:

```r
min_N <- 150
```

Set `min_N <- NA` to disable sample-size filtering.

#### Minimum marker pairs per bin

```r
min_pairs_per_bin <- 1
```

Increasing this value can remove unstable bins with very few marker pairs.

#### Chromosome name conversion

```r
convert_chr_names <- TRUE
```

This converts TASSEL-style chromosome labels such as `A`, `B`, `C`, ..., `J` to `Chr1`, `Chr2`, `Chr3`, ..., `Chr10`.

### Raw point cloud settings

These settings apply only to `GWAS_LD_scatter_graph.R`.

```r
show_raw_point_cloud <- TRUE
max_raw_points_to_plot <- Inf
raw_point_seed <- 123
raw_point_alpha <- 0.40
raw_point_size <- 0.30
raw_point_color <- "grey45"
```

To plot every raw pairwise LD value:

```r
max_raw_points_to_plot <- Inf
```

To reduce output file size and plotting time, use a fixed subset:

```r
max_raw_points_to_plot <- 300000
```

The subset is reproducible because the script uses:

```r
raw_point_seed <- 123
```

## Output file descriptions

### GWAS outputs

#### `<dataset>_manhattan.png` and `<dataset>_manhattan.svg`

Publication-ready Manhattan plots.

The SVG version is useful for final editing in vector-graphics software such as Inkscape, Adobe Illustrator, or Affinity Designer.

#### `<dataset>_qq.png` and `<dataset>_qq.svg`

QQ plots with confidence envelopes and genomic inflation factor.

#### `<dataset>_cleaned_chromosome_markers.csv`

Filtered GWAS marker table used for plotting and summary statistics.

#### `<dataset>_top_hits.csv`

Top markers sorted by ascending p-value.

#### `<dataset>_chromosome_marker_summary.csv`

Per-chromosome marker summary containing:

```text
Chr
n_markers
min_pos
max_pos
min_p
max_log10p
```

#### `CTI_summary_statistics.csv`

Overall summary table containing, for each dataset:

```text
dataset
input_file
n_markers
n_chr
min_p
max_log10p
lambda_gc
bonferroni_p
suggestive_p
bonferroni_log10p
suggestive_log10p
qq_confidence_level
top_marker
top_chr
top_pos
top_p
```

If the scripts are used for traits other than CTI, this output filename can be renamed in the script.

### LD outputs

#### `LD_decay_binned_values.csv`

Binned LD values used for smoothed or fitted LD decay curves.

Columns include:

```text
chr
dist_kb
mean_r2
median_r2
n_pairs
```

#### `LD_decay_<method>.png` and `LD_decay_<method>.svg`

LD decay plots for the selected method.

For `GWAS_LD_scatter_graph.R`, the plot also includes the raw grey pairwise LD point cloud.

#### `LD_decay_<method>_curve_values.csv`

Fitted or smoothed LD curve values used for plotting.

#### `LD_decay_<method>_table.csv`

Final LD decay table by chromosome, including decay distances at the selected r² thresholds.

For the `hill_weir` method, the table may also include:

```text
C
fit_status
LD_half_decay_kb
half_decay_threshold
```

A `Mean` row is added when enough chromosomes have valid decay estimates.

#### `LD_decay_<method>_decay_diagnostics.csv`

Diagnostic table reporting whether a decay threshold was successfully crossed for each chromosome and threshold.

Possible status values include:

```text
ok
ok_no_interpolation
too_few_points
below_first_point
always_below_threshold
not_reached_within_max_dist
no_downward_crossing
```

## LD decay methods

### Raw bin mean

```r
decay_method <- "raw_bin_mean"
```

Uses the mean r² value within each distance bin without additional smoothing.

This is transparent but can be noisy.

### Rolling mean

```r
decay_method <- "rolling_mean"
```

Calculates the binned mean r² and smooths it with a moving average.

Relevant setting:

```r
rolling_window <- 7
```

### LOESS

```r
decay_method <- "loess"
```

Fits a LOESS curve to binned mean r² values.

Relevant settings:

```r
loess_span <- 0.25
loess_degree <- 2
```

Smaller `loess_span` values follow local fluctuations more closely. Larger values produce smoother curves.

### Hill-Weir / Remington model

```r
decay_method <- "hill_weir"
```

Fits the expected r² decay model described by Hill and Weir and commonly used in LD decay analyses following Remington et al.

This method requires the `N` column in the TASSEL LD export.

Relevant settings:

```r
hill_weir_max_points_per_chr <- 100000
hill_weir_seed <- 123
hill_weir_start_values <- c(1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 0.01, 0.1)
hill_weir_prediction_step_kb <- 1
```

For large LD files, the script can randomly sample marker pairs per chromosome for model fitting while keeping the sampling reproducible.

## Recommended workflow

1. Export GWAS results from TASSEL.
2. Check that the GWAS table contains at least `Trait`, `Marker`, `Chr`, `Pos`, and `p`.
3. Edit `input_files` in `GWAS_graph.R`.
4. Run:

```bash
Rscript GWAS_graph.R
```

5. Inspect the cleaned marker tables and summary statistics.
6. Export LD results from TASSEL.
7. Check that the LD table contains `Locus1`, `Locus2`, `Dist_bp`, and `R^2`.
8. Edit `input_csv` in either LD script.
9. Select the LD decay method and thresholds.
10. Run either:

```bash
Rscript GWAS_LD_graph.R
```

or:

```bash
Rscript GWAS_LD_scatter_graph.R
```

11. Use the SVG outputs for final publication figure editing.

## Notes and limitations

* The scripts currently use file paths defined directly inside the R scripts. They do not yet use command-line arguments.
* TASSEL chromosome labels are expected to be compatible with the chromosome ordering defined in the script.
* For GWAS plotting, only chromosomes `A` to `J` are retained by default.
* For LD plotting, only intra-chromosomal marker pairs are used.
* The Hill-Weir model requires valid sample-size information in the `N` column.
* Plotting all raw LD points can generate very large SVG files. For large datasets, use a finite value for `max_raw_points_to_plot`.
* LD decay distances are reported as `NA` if the selected threshold is not crossed within the analyzed distance range or if the curve is already below the threshold at the first evaluable point.
* The scripts are intended for post-processing and visualization of TASSEL outputs, not for performing GWAS or LD calculations themselves.

## Example

Example command-line execution:

```bash
Rscript GWAS_graph.R
Rscript GWAS_LD_scatter_graph.R
```

Example LD settings for a publication-style LD decay plot with raw point cloud:

```r
max_dist_kb <- 2000
bin_size_kb <- 1
thresholds <- c(0.2, 0.1)
decay_method <- "hill_weir"
show_raw_point_cloud <- TRUE
max_raw_points_to_plot <- Inf
raw_point_alpha <- 0.40
raw_point_size <- 0.30
raw_point_color <- "grey45"
```

## Citation

If you use these scripts, please cite the associated publication:

```text
[Add manuscript citation here after publication]
```

The LD decay model implemented in the Hill-Weir option follows:

```text
Hill WG, Weir BS. 1988. Variances and covariances of squared linkage disequilibria in finite populations. Theoretical Population Biology.

Remington DL et al. 2001. Structure of linkage disequilibrium and phenotypic associations in the maize genome. PNAS 98:11479–11484.
```

## License

GPL-3.0

## Contact

For questions about the scripts, input formatting, or reproducibility of the maize GWAS figures, please contact the corresponding author of the associated manuscript.
