# Termite-evo-devo

Key code developed for **"Tracking major evolutionary transitions to superorganisms"**, including quantification of between-caste divergence, evaluation of caste gene-regulatory similarity, measurement of between-caste distance across developmental stages, and principal-component projection of cross-species caste developmental trajectories.

## Repository contents

| File | Description |
|---|---|
| `Supplementary_codes.R` | R functions used in the analyses (see [Functions](#functions)) |
| `run_demo.R` | Script that runs all analyses on the demo dataset |
| `termite_exp_demo.Rdata` | Demo expression data: head–prothorax RNA-seq of four termite species |
| `termite_gene_info_demo.Rdata` | Demo gene annotation and orthogroups for the four species |
| `LICENSE` | MIT License |

---

## 1. System requirements

### Software dependencies

| Software | Version tested | Source |
|---|---|---|
| R | 4.3.3 | https://www.r-project.org |
| dplyr | 1.1.4 | CRAN |
| EDASeq | 2.36.0 | Bioconductor (needed only for `combine_exp_data()`) |
| sva | 3.50.0 | Bioconductor (needed only for `reverse_pca_multi()`) |
| ggplot2 | 3.4.4 | CRAN (needed only for `plot_single_candidate_dmel()`) |

All other functions used (`lm`, `svd`, `dist`, `cor`, base graphics) are part of base R.

### Operating systems tested

- Ubuntu 24.04 LTS (x86_64), R 4.3.3
- macOS 27 (Apple silicon), R 4.6.1 -->

The code is plain R and is expected to run on any operating system supported by R ≥ 4.0 (Linux, macOS, Windows).

### Hardware

No non-standard hardware is required. The demo runs on a standard desktop or laptop computer with ≥ 4 GB RAM.

---

## 2. Installation guide

### Instructions

1. Install R (≥ 4.0) from https://cran.r-project.org.
2. Install the required packages from within R:

   ```r
   install.packages(c("dplyr", "ggplot2"))

   if (!requireNamespace("BiocManager", quietly = TRUE))
       install.packages("BiocManager")
   BiocManager::install(c("EDASeq", "sva"))
   ```

3. Download this repository:

   ```bash
   git clone https://github.com/BitaoQiu/Termite-evo-devo.git
   cd Termite-evo-devo
   ```

   or download the ZIP file from GitHub (**Code → Download ZIP**) and unzip it.

No compilation is needed; the functions are loaded with `source("Supplementary_codes.R")`.

### Typical install time

About 5–10 minutes on a normal desktop computer for the R packages (most of this is installing `EDASeq`, `sva` and their Bioconductor dependencies). Installing R itself, if needed, adds a few minutes.

---

## 3. Demo

### Demo data

The demo data are head–prothorax RNA-seq data for four termite species, in two files.

**`termite_exp_demo.Rdata`** contains `termite_exp_data_zscore_HP`, a named list with one entry per species:

| Species | Code | Genes | Samples | Castes / stages |
|---|---|---|---|---|
| *Zootermopsis nevadensis* | `Znev` | 15,845 | 25 | workers (larva, early, middle), nymphs (middle, late) |
| *Mastotermes darwiniensis* | `Mdar` | 16,550 | 28 | workers (larva, middle), nymphs (middle, late, final) |
| *Reticulitermes grassei* | `Rgra` | 18,103 | 46 | workers (early, middle), nymphs (middle, late), adult kings and queens |
| *Macrotermes barneyi* | `Mbar` | 19,703 | 51 | small and large workers, small and large soldiers (middle to adult), nymphs (early to adult) |

Each species entry contains:

- `abundance` — normalized expression matrix, log2(TPM) scale (genes × samples)
- `sampleInfo` — sample metadata, one row per sample in the same order as the matrix columns; key columns are `Sample_ID`, `Caste`, `Stage`, `Stage_Caste`, `BatchDate`, and `Species`

**`termite_gene_info_demo.Rdata`** contains:

- `gene_info_demo` — gene annotation per species: `GeneID`, `Description`, `Dmel_name` and `Dmel_bitScore` (*Drosophila melanogaster* best hit), GO/KEGG/Reactome terms, and gene age
- `ortho_info_demo` — orthogroups in OrthoFinder format (13,832 orthogroups): an `Orthogroup` column plus one column per species listing that species' genes in the group, comma-separated; 8,772 orthogroups are one-to-one across the four species

### Instructions to run on data

From the repository folder, run:

```bash
Rscript run_demo.R
```

or, inside an R session started in the repository folder:

```r
source("run_demo.R")
```

The script runs:

0. **Cross-species data set** — one-to-one orthologs are extracted from the orthogroups, and the four species are combined into one ortholog-matched data set.
1. **Between-caste divergence** — per-gene Hedges' g for workers vs. nymphs, adjusted for sequencing batch, in each species; summarized as the median |g| across genes. Worker castes (e.g. small and large workers in `Mbar`) are pooled.
2. **Similarity of caste-biased expression between species** — angular distance between the species' worker-vs-nymph Hedges' g profiles over one-to-one orthologs, and a list of orthologs biased in the same direction in all species.
3. **Between-caste distance across development** — Manhattan distance (per gene) of each nymph sample to middle-stage workers, by nymph stage, with the within-worker distance as a baseline.
4. **Cross-species PCA projection** — PC space built from *Z. nevadensis*; the other three species are projected onto it after ComBat alignment (species as batch, reference unchanged).
5. **Illustration** — expression of *Kr-h1* (Krüppel homolog 1) across stages and castes in each species.

### Expected output

The script prints progress and summary tables to the console and writes the following to `demo_output/`:

| File | Content |
|---|---|
| `0_one2one_orthologs.csv` | One-to-one ortholog table used for cross-species analyses |
| `1_hedges_g_<species>.csv` | Per-gene coefficient, residual SD, Cohen's d, Hedges' g, and gene annotation |
| `1_hedges_g_summary.csv` | Number of genes and median \|g\| per species |
| `2_angular_distance_species.csv` | Angular distance between species' caste-bias profiles (0 = same direction, 0.5 = orthogonal, 1 = opposite) |
| `2_consistent_caste_bias.csv` | Orthologs biased in the same direction in all species, with annotation |
| `3_sample_distance_<species>.csv` | Distance of every nymph sample to every middle-stage worker sample, with metadata |
| `3_sample_distance_summary.csv` | Median distance per species and nymph stage, and the within-worker baseline |
| `4_reverse_pca.csv` | PC1–PC5 scores for reference (`Znev`) and projected samples |
| `4_reverse_pca_group_means.csv` | Mean PC1/PC2 per species, stage and caste |
| `5_Kr-h1_expression.pdf` | *Kr-h1* expression by stage and caste, one page per species |
| `demo_figures.pdf` | Summary figure: (A) Hedges' g distributions, (B) angular distance between species, (C) nymph-to-worker distances by stage, (D) cross-species PCA projection |

Expected console output (values may differ in the last decimal places across platforms; ComBat progress messages omitted):

```
[0] One-to-one orthologs and combined data set
8772 one-to-one ortholog groups across Znev, Mdar, Rgra, Mbar.
Combined: 8772 ortholog groups x 150 samples

[1] Hedges' g (Worker vs. Nymph, batch-adjusted)
 Species n_genes median_abs_g
    Znev   10360    0.3158563
    Mdar    9383    0.5257383
    Rgra   10944    0.6259475
    Mbar   10924    0.6480816

[2] Angular distance between species' caste-biased expression profiles
Orthologs with g in all species: 6021
      Znev  Mdar  Rgra  Mbar
Znev 0.000 0.500 0.481 0.493
Mdar 0.500 0.000 0.435 0.499
Rgra 0.481 0.435 0.000 0.400
Mbar 0.493 0.499 0.400 0.000

Orthologs biased in the same direction in all species: 1025; top 10:
 Orthogroup  Znev  Mdar  Rgra  Mbar   Bias Dmel_name                              Description
  OG0007340  1.29  3.67  3.51  1.49 Worker     Kr-h1                                     <NA>
  OG0008488  1.21  1.99  2.31  2.06 Worker      <NA>                                     <NA>
  OG0003192 -1.37 -1.09 -1.15 -3.36  Nymph      Tps1                                     <NA>
  OG0007313 -1.06 -1.35 -6.08 -2.81  Nymph      <NA>                                     <NA>
  OG0009556 -1.30 -1.11 -1.61 -1.00  Nymph     Ntan1              N-terminal amidohydrolase 1
  OG0006702  1.11  0.88  1.43  1.27 Worker     VAChT vesicular acetylcholine transporter-like
  OG0009166 -0.86 -5.93 -4.03 -3.00  Nymph      Arr2                               arrestin 2
  OG0000910 -1.18 -2.16 -2.73 -0.83  Nymph      Snoo                             Sno oncogene
  OG0007095 -0.98 -1.58 -1.93 -0.82  Nymph  Oscillin glucosamine-6-phosphate isomerase Oscill
  OG0008526  0.88  0.80  1.74  1.24 Worker      Syt7                    synaptotagmin-15-like

[3] Distance of nymph samples to middle-stage workers
 Species  Stage median_distance within_worker_baseline
    Znev Middle       0.3277889              0.3228932
    Znev   Late       0.3973463              0.3228932
    Mdar Middle       0.4095333              0.3406655
    Mdar   Late       0.4362442              0.3406655
    Mdar  Final       0.6929000              0.3406655
    Rgra Middle       0.3804718              0.2849457
    Rgra   Late       0.3593370              0.2849457
    Mbar  Early       0.5216510              0.3605829
    Mbar Middle       0.3843661              0.3605829
    Mbar   Late       0.5313572              0.3605829
    Mbar  Adult       0.7877894              0.3605829

[4] Reference-based PCA projection (reference: Znev)
Projected 125 samples onto the PC space of 25 Znev samples.

[5] Expression plots of Kr-h1
Znev: plotted gene-Znev_012015
Mdar: plotted gene-Mdar_013062
Rgra: plotted gene-Rgra_003466
Mbar: plotted gene-Mbar_008534

Done. Output written to ./demo_output/ (run time: ~130 s)
```

### Expected run time

About 2 minutes on a normal desktop computer (measured at 129 s on a 2-core machine). Most of the time is spent fitting the per-gene linear models for Hedges' g.

---

## 4. Instructions for use

Load the functions in R:

```r
source("Supplementary_codes.R")
```

### Input format

All functions take an expression matrix with **genes as rows and samples as columns** (column names = sample IDs), plus a sample-metadata data frame with one row per sample. Metadata rows must be in the same order as the matrix columns for `cal_hedges_g()`; the other functions match samples by ID. Expression values should be normalized and on a log-like scale (e.g. log-transformed TPM or variance-stabilized counts).

### Functions

#### `get_one2one_orthologs()` — one-to-one orthologs

Extracts orthogroups with exactly one gene in every species from an OrthoFinder-style table (e.g. `Orthogroups.tsv`):

```r
ortho   <- read.delim("Orthogroups.tsv", check.names = FALSE)  # columns: Orthogroup, <species>...
one2one <- get_one2one_orthologs(ortho, species_list = c("Znev", "Mdar", "Rgra"))
```

Species column names must match the species codes used in the expression data. Groups that share a gene with another retained group are removed with a warning.

#### `combine_exp_data()` — ortholog-matched multi-species data set

Combines per-species expression data into one data set with rows matched by one-to-one orthologs, as input for cross-species analyses:

```r
# exp_list: named list per species, each with abundance and sampleInfo (optionally counts)
# one2one:  output of get_one2one_orthologs()
combined <- combine_exp_data(exp_data        = exp_list,
                             ortho_info      = one2one,
                             quantile.method = "upper",   # upper-quartile normalization
                             species_list    = c("Znev", "Mdar", "Rgra"))
```

The output holds `abundance` (ortholog groups present in all species), `abundance.quantile` (between-sample normalized across species), `counts` (if given for all species), and the merged `sampleInfo`. Each species' `sampleInfo` needs the columns `Sample_ID`, `Ind`, `Caste`, `Colony`, `Stage`, `Stage_Caste`, `Tissue`, `Species`, `BatchDate` and `Batch`, in the same order as the matrix columns.

#### `cal_hedges_g()` — between-caste divergence

Hedges' g for one gene from a linear model that can include covariates such as batch. Apply it to all genes and summarize with the median:

```r
res <- t(sapply(rownames(abundance), cal_hedges_g,
                abundance  = abundance,
                sampleInfo = sampleInfo,
                model      = "exp ~ Caste + BatchDate",
                co_var     = "CasteWorker",          # coefficient name as in coef(fit)
                ref_levels = list(Caste = "Nymph")))
median(abs(res[, "Hedges_g"]), na.rm = TRUE)
```

The response variable in `model` must be named `exp`. A gene returns `NA` when the contrast is not estimable (for example, when caste is fully confounded with batch).

#### `angular_distance_matrix_cols()` — profile similarity

Angular distance between all pairs of columns of a matrix, for example log fold-change profiles of caste-biased expression from several species:

```r
# lfc: genes (orthologs) x contrasts, e.g. one column per species
d <- angular_distance_matrix_cols(lfc)
```

For cross-species comparisons, rows must be ortholog-matched genes in the same order for all columns. Missing values are treated as 0 (no change).

#### `cal_sample_dist()` — between-caste distance across stages

Distance of each reference sample (e.g. nymphs at each stage) to a set of target samples (e.g. middle-stage workers):

```r
df <- cal_sample_dist(exp_data     = list(abundance = abundance, sampleInfo = sampleInfo),
                      target_caste = "Worker",
                      target_stage = "Middle",
                      ref_caste    = "Nymph",
                      ref_stages   = c("Early", "Middle", "Late", "Adult"),
                      dist_fun     = "man")   # or "s", "p", "k", "angular"
```

`sampleInfo` needs the columns `Sample_ID`, `Caste`, `Stage` and `Species` (and `BatchDate` if `subset_BatchDate` is used). Only genes with median abundance > 1 are used.

#### `reverse_pca_multi()` — reference-based PCA projection

Builds the PC space from a reference group (e.g. one species) and projects the other groups onto it:

```r
# combined: output of combine_exp_data()
info <- combined$sampleInfo
rownames(info) <- info$Sample_ID
groups <- split(info$Sample_ID, info$Species)   # list(Znev = ..., Mdar = ..., Rgra = ...)
pca <- reverse_pca_multi(exp_data      = combined$abundance,
                         sample.info   = info,
                         sample_groups = groups,
                         ref_batch     = "Znev",
                         batch_col     = "Species",
                         correct_batch = TRUE,    # ComBat alignment
                         combat_ref    = TRUE,    # keep the reference unchanged
                         filter        = TRUE,    # top-variance genes in the reference
                         n_gene        = 512,
                         n_pc          = 5)
```

The output contains PC scores for all samples joined with their metadata; the column `set_` marks reference (`ref`) and projected (`target`) samples. `run_demo.R` contains complete worked examples of all functions.

#### `plot_single_candidate_dmel()` — expression of one gene (illustration)

Plots the expression of the gene whose *Drosophila* best hit is `dmel_name` (the hit with the highest bit score if several genes share the name), as box plots or points:

```r
p <- plot_single_candidate_dmel("Kr-h1",
                                exp_data  = exp_list$Znev,          # abundance + sampleInfo
                                gene_info = gene_info_demo$Znev,    # GeneID, Dmel_name, Dmel_bitScore, Description
                                x = "Stage_Caste", col = "Caste",
                                plot_type = "boxplot")              # or "point"
print(p)   # a ggplot object; further ggplot2 layers can be added
```

---
## License

This project is licensed under the MIT License; see `LICENSE`.
