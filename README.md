# Termite-evo-devo

Key code developed for **"Tracking major evolutionary transitions to superorganisms"**, including quantification of between-caste divergence, evaluation of caste gene-regulatory similarity, measurement of between-caste distance across developmental stages, and principal-component projection of cross-species caste developmental trajectories.

## Repository contents

| File | Description |
|---|---|
| `termite_evo_devo_functions.R` | Key R functions developed in the analyses (see [Functions](#functions)) |
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
| limma | 3.58.1 | Bioconductor (used in the demo for caste DEGs; installed with `sva`) |
| ggplot2 | 3.4.4 | CRAN (needed only for `plot_single_candidate_dmel()`) |

All other functions used (`lm`, `svd`, `dist`, `cor`, base graphics) are part of base R.

### Operating systems tested

- Ubuntu 24.04 LTS (x86_64), R 4.3.3
- macOS 27 (Apple silicon), R 4.6.1

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
   BiocManager::install(c("EDASeq", "sva", "limma"))
   ```

3. Download this repository:

   ```bash
   git clone https://github.com/BitaoQiu/Termite-evo-devo.git
   cd Termite-evo-devo
   ```

   or download the ZIP file from GitHub (**Code → Download ZIP**) and unzip it.

No compilation is needed; the functions are loaded with `source("termite_evo_devo_functions.R")`.

### Typical install time

About 5–10 minutes on a normal desktop computer for the R packages (most of this is installing `EDASeq`, `sva` and their Bioconductor dependencies). Installing R itself, if needed, adds a few minutes.

---

## 3. Demo

### Demo data

The demo data are head–prothorax RNA-seq data for four termite species, in two files.

**`termite_exp_demo.Rdata`** contains `termite_exp_data_HP`, a named list with one entry per species:

| Species | Code | Genes | Samples | Castes / stages |
|---|---|---|---|---|
| *Cryptotermes secundus* | `Csec` | 20,736 | 27 | workers (larva, early, middle), nymphs (middle, late) |
| *Mastotermes darwiniensis* | `Mdar` | 16,550 | 28 | workers (larva, middle), nymphs (early, late, final) |
| *Reticulitermes grassei* | `Rgra` | 18,103 | 46 | workers (early, middle), nymphs (early, late), adult kings and queens |
| *Macrotermes barneyi* | `Mbar` | 19,703 | 51 | small and large workers, small and large soldiers (middle to adult), nymphs (early to adult) |

Each species entry contains:

- `abundance` — normalized expression matrix (batch corrected), log2(TPM) scale (genes × samples)
- `sampleInfo` — sample metadata, one row per sample in the same order as the matrix columns; columns are `Sample_ID`, `Ind`, `Caste`, `Stage`, `Stage_Caste`, `Tissue`, `Colony` and `BatchDate` (the demo adds `Species` from the list name)

**`termite_gene_info_demo.Rdata`** contains:

- `gene_info_demo` — gene annotation per species: `GeneID`, `Description`, `Dmel_name` and `Dmel_bitScore` (*Drosophila melanogaster* best hit), GO/KEGG/Reactome terms, and gene age
- `ortho_info_demo` — orthogroups in OrthoFinder format (13,832 orthogroups): an `Orthogroup` column plus one column per species listing that species' genes in the group, comma-separated; 8,666 orthogroups are one-to-one across the four species

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
1. **Between-caste divergence** — per-gene Hedges' g for workers vs. nymphs in each species; summarized as the median |g| across genes. Positive g means higher expression in workers. Each species' contrast uses these `Stage_Caste` groups:

   | Species | Workers | Nymphs |
   |---|---|---|
   | `Csec` | Middle_Worker | Late_Nymph |
   | `Mdar` | Middle_Worker | Late_Nymph |
   | `Rgra` | Middle_Worker | Late_Nymph |
   | `Mbar` | Adult_Small_worker | Late_Nymph |

   Because the demo abundance data are already batch corrected, the model contains caste group only (`exp ~ Group`). For data that are not batch corrected, batch can be added as a covariate (see [`cal_hedges_g()`](#cal_hedges_g--between-caste-divergence)).

2. **Similarity of caste-biased expression between species** — angular distance between the species' worker-vs-nymph log2 fold-change (LFC) profiles, over one-to-one orthologs that are caste DEGs in at least one species; |LFC| is capped at 3 so that a few extreme genes do not dominate. Also lists orthologs biased in the same direction in all species. The LFC contrasts differ from step 1:

   | Species | Workers | Nymphs |
   |---|---|---|
   | `Csec` | Middle_Worker | Middle_Nymph + Late_Nymph |
   | `Mdar` | Middle_Worker | Early_Nymph |
   | `Rgra` | Middle_Worker | Early_Nymph |
   | `Mbar` | Middle_Worker | Early_Nymph |

   *This step is a simplified illustration:* LFC is the difference in mean log2(TPM) between the two groups, and DEGs are called with limma's moderated t-test (BH-adjusted p < 0.05) on the demo data. In the manuscript, LFCs and DEGs come from formal differential expression analysis; with real data, those LFCs are passed to `angular_distance_matrix_cols()` in the same way.
3. **Between-caste distance across development** — Manhattan distance (per gene) of each nymph sample to middle-stage workers, by nymph stage, with the within-worker distance as a baseline.
4. **Cross-species PCA projection on caste genes** — the top 200 caste DEGs in *C. secundus* (all nymphs vs. Middle_Worker, Welch t-test on abundance, among one-to-one orthologs) are selected; the PC space is built from *C. secundus* (the reference species in the manuscript) using only these genes, and the other three species are projected onto it after ComBat alignment (species as batch, reference unchanged). In all species, only Early_Nymph, Middle_Nymph and Late_Nymph (whichever are present), Middle_Worker and Middle_Small_worker samples are used, plus Adult_Worker and Adult_Small_worker in `Mbar`.
5. **Illustration** — expression of *Kr-h1* (Krüppel homolog 1) and *Eip93F* (Ecdysone-induced protein 93F) across stages and castes in each species, with `plot_single_candidate_dmel()`.

### Expected output

The script prints progress and summary tables to the console and writes the following to `demo_output/`:

| File | Content |
|---|---|
| `1_hedges_g_summary.csv` | Contrast, sample sizes, number of genes and median \|g\| per species |
| `2_angular_distance_species.csv` | Angular distance between species' caste-bias profiles (0 = same direction, 0.5 = orthogonal, 1 = opposite) |
| `2_angular_distance_tree.pdf` | Hierarchical clustering of species by angular distance (`hclust`) |
| `2_consistent_caste_bias.csv` | Orthologs biased in the same direction in all species, with capped LFC per species, number of species in which the gene is a DEG, and annotation |
| `3_sample_distance_summary.csv` | Median distance per species and nymph stage, and the within-worker baseline |
| `4_csec_caste_DEGs.csv` | The 200 *C. secundus* caste DEGs used for the PCA: orthogroup, mean difference (worker − nymph), t, p, bias, gene ID and *Drosophila* name |
| `4_reverse_pca.pdf` | PC1 vs. PC2 per species, coloured by caste and shaped by stage |
| `5_Kr-h1_expression.pdf`, `5_Eip93F_expression.pdf` | *Kr-h1* and *Eip93F* expression by stage and caste, one page per species |

Expected console output (values may differ in the last decimal places across platforms; ComBat progress messages omitted):

```
[0] One-to-one orthologs and combined data set
8666 one-to-one ortholog groups across Csec, Mdar, Rgra, Mbar.
Combined: 8666 ortholog groups x 152 samples

[1] Hedges' g (Worker vs. Nymph)
 Species                          Contrast n_workers n_nymphs n_genes median_abs_g
    Csec      Middle_Worker vs. Late_Nymph        10        3   11458    0.6743113
    Mdar      Middle_Worker vs. Late_Nymph        11        3    9366    0.8555363
    Rgra      Middle_Worker vs. Late_Nymph        15        4   10931    0.7010929
    Mbar Adult_Small_worker vs. Late_Nymph         4        4   10974    1.4543576

[2] Angular distance between species' caste-biased expression profiles
Caste DEGs per species (one-to-one orthologs):
 Species                                  Contrast n_DEG
    Csec Middle_Worker vs. Middle_Nymph+Late_Nymph   104
    Mdar             Middle_Worker vs. Early_Nymph  1793
    Rgra             Middle_Worker vs. Early_Nymph  1215
    Mbar             Middle_Worker vs. Early_Nymph   778
Orthologs that are DEGs in at least one species: 3239
      Csec  Mdar  Rgra  Mbar
Csec 0.000 0.483 0.424 0.513
Mdar 0.483 0.000 0.424 0.470
Rgra 0.424 0.424 0.000 0.497
Mbar 0.513 0.470 0.497 0.000

Orthologs biased in the same direction in all species: 680; top 10:
 Orthogroup  Csec  Mdar  Rgra  Mbar   Bias n_species_DEG Dmel_name                              Description
  OG0003120 -3.00 -3.00 -3.00 -3.00  Nymph             4      Lsp2                           hexamerin-like
  OG0007313 -3.00 -3.00 -3.00 -3.00  Nymph             4      <NA>                                     <NA>
  OG0009994 -3.00 -3.00 -3.00 -2.85  Nymph             2      Lsp2                           hexamerin-like
  OG0007185 -1.32 -1.35 -1.94 -3.00  Nymph             4        cd             peroxidasin homolog cardinal
  OG0009166 -1.32 -3.00 -1.41 -2.07  Nymph             4      Arr2                               arrestin 2
  OG0007811 -1.51 -2.09 -1.08 -2.27  Nymph             3        nw                                   narrow
  OG0007340  1.02  1.16  1.13  1.19 Worker             3     Kr-h1                                     <NA>
  OG0000879 -0.92 -1.11 -1.12 -1.01  Nymph             3       Jhe           juvenile hormone esterase-like
  OG0002392 -3.00 -3.00 -0.87 -3.00  Nymph             3      <NA>                                     <NA>
  OG0005137 -0.74 -1.11 -0.78 -1.11  Nymph             3       chp leucine rich repeat containing G protein

[3] Distance of nymph samples to middle-stage workers
 Species  Stage median_distance within_worker_baseline
    Csec Middle       0.4004395              0.3057668
    Csec   Late       0.3478501              0.3057668
    Mdar  Early       0.4095333              0.3406655
    Mdar   Late       0.4362442              0.3406655
    Mdar  Final       0.6929000              0.3406655
    Rgra  Early       0.3804718              0.2849457
    Rgra   Late       0.3593370              0.2849457
    Mbar  Early       0.5216510              0.3605829
    Mbar Middle       0.3843661              0.3605829
    Mbar   Late       0.5313572              0.3605829
    Mbar  Adult       0.7877894              0.3605829

[4] Reference-based PCA projection on Csec caste DEGs (reference: Csec)
Top 200 Csec caste DEGs (Middle_Worker, n = 10 vs. Late_Nymph+Middle_Nymph, n = 7): 93 worker-biased, 107 nymph-biased; max p = 0.0042
      
       Early_Nymph Middle_Nymph Late_Nymph Middle_Worker Middle_Small_worker Adult_Worker Adult_Small_worker
  Csec           0            4          3            10                   0            0                  0
  Mdar           6            0          3            11                   0            0                  0
  Rgra           4            0          4            15                   0            0                  0
  Mbar           4            4          4             4                   4            4                  4
Projected 71 samples onto the PC space of 17 Csec samples.

[5] Expression plots of Kr-h1 and Eip93F 
Kr-h1 Csec: plotted gene-Csec_013894
Kr-h1 Mdar: plotted gene-Mdar_013062
Kr-h1 Rgra: plotted gene-Rgra_003466
Kr-h1 Mbar: plotted gene-Mbar_008534
Eip93F Csec: plotted gene-Csec_019594
Eip93F Mdar: plotted gene-Mdar_002584
Eip93F Rgra: plotted gene-Rgra_011182
Eip93F Mbar: plotted gene-Mbar_011384

Done. Output written to ./demo_output/ (run time: ~110 s)
```

### Expected run time

About 2 minutes on a normal desktop computer (measured at 100–140 s on a 2-core machine). Most of the time is spent fitting the per-gene linear models for Hedges' g.

---

## 4. Instructions for use

Load the functions in R:

```r
source("termite_evo_devo_functions.R")
```

### Input format

All functions take an expression matrix with **genes as rows and samples as columns** (column names = sample IDs), plus a sample-metadata data frame with one row per sample. Metadata rows must be in the same order as the matrix columns for `cal_hedges_g()`; the other functions match samples by ID. Expression values should be normalized and on a log-like scale (e.g. log-transformed TPM or variance-stabilized counts).

### Functions

#### `get_one2one_orthologs()` — one-to-one orthologs

Extracts orthogroups with exactly one gene in every species from an OrthoFinder-style table (e.g. `Orthogroups.tsv`):

```r
ortho   <- read.delim("Orthogroups.tsv", check.names = FALSE)  # columns: Orthogroup, <species>...
one2one <- get_one2one_orthologs(ortho, species_list = c("Csec", "Mdar", "Rgra"))
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
                             species_list    = c("Csec", "Mdar", "Rgra"))
```

The output holds `abundance` (ortholog groups present in all species), `abundance.quantile` (between-sample normalized across species), `counts` (if given for all species), and the merged `sampleInfo`. Each species' `sampleInfo` needs the columns `Sample_ID`, `Caste`, `Stage` and `Stage_Caste`, in the same order as the matrix columns. `Species` is taken from the list name if absent; `Ind`, `Colony`, `Tissue`, `BatchDate` and `Batch` are kept if present.

#### `cal_hedges_g()` — between-caste divergence

Hedges' g for one gene from a linear model that can include covariates such as batch (omit them if the data are already batch corrected, as in the demo). Apply it to all genes and summarize with the median:

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
# combined:  output of combine_exp_data()
# caste_degs: orthogroup IDs of the genes to use (e.g. top caste DEGs in the
#             reference species); omit `genes` to use all genes
info <- combined$sampleInfo
rownames(info) <- info$Sample_ID
groups <- split(info$Sample_ID, info$Species)   # list(Csec = ..., Mdar = ..., Rgra = ...)
pca <- reverse_pca_multi(exp_data      = combined$abundance,
                         sample.info   = info,
                         sample_groups = groups,
                         ref_batch     = "Csec",
                         genes         = caste_degs,
                         batch_col     = "Species",
                         correct_batch = TRUE,    # ComBat alignment
                         combat_ref    = TRUE,    # keep the reference unchanged
                         filter        = FALSE,   # TRUE: top-variance genes in the reference instead
                         n_pc          = 5)
```

The output contains PC scores for all samples joined with their metadata; the column `set_` marks reference (`ref`) and projected (`target`) samples. `run_demo.R` contains complete worked examples of all functions.

#### `plot_single_candidate_dmel()` — expression of one gene (illustration)

Plots the expression of the gene whose *Drosophila* best hit is `dmel_name` (the hit with the highest bit score if several genes share the name), as box plots or points:

```r
p <- plot_single_candidate_dmel("Kr-h1",
                                exp_data  = exp_list$Csec,          # abundance + sampleInfo
                                gene_info = gene_info_demo$Csec,    # GeneID, Dmel_name, Dmel_bitScore, Description
                                x = "Stage_Caste", col = "Caste",
                                plot_type = "boxplot")              # or "point"
print(p)   # a ggplot object; further ggplot2 layers can be added
```

---
## License

This project is licensed under the MIT License; see `LICENSE`.
