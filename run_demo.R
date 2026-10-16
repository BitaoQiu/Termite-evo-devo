# =============================================================================
# Demo: runs the analyses in termite_evo_devo_functions.R on the bundled demo data:
#   termite_exp_demo.Rdata        head-prothorax RNA-seq of four termite species
#   termite_gene_info_demo.Rdata  gene annotation and orthogroups
#
# Usage (from the repository folder):
#   Rscript run_demo.R
#
# Output: tables and figures are written to ./demo_output/
# =============================================================================

t_start <- Sys.time()
set.seed(1)
options(width = 120)

source("termite_evo_devo_functions.R")   # loads dplyr and defines the functions

out_dir <- "demo_output"
dir.create(out_dir, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Load demo data
#   termite_exp_demo.Rdata
#     termite_exp_data_HP : named list, one element per species, each
#                                  holding abundance (genes x samples) and
#                                  sampleInfo (one row per sample, same order)
#   termite_gene_info_demo.Rdata
#     gene_info_demo             : gene annotation per species (GeneID, Description,
#                                  Dmel_name = Drosophila best hit, Dmel_bitScore, ...)
#     ortho_info_demo            : orthogroups (OrthoFinder format)
# -----------------------------------------------------------------------------
load("termite_exp_demo.Rdata")
load("termite_gene_info_demo.Rdata")
exp_list <- termite_exp_data_HP
species  <- c("Csec", "Mdar", "Rgra", "Mbar")
# Species label for each sample (used when combining species and in summaries)
for (sp in species) exp_list[[sp]]$sampleInfo$Species <- sp
for (sp in species) {
  cat(sprintf("%s: %d genes x %d samples\n", sp,
              nrow(exp_list[[sp]]$abundance), ncol(exp_list[[sp]]$abundance)))
}

# Worker castes differ between species (e.g. small and large workers in Mbar)
worker_castes <- c("Worker", "Small_worker", "Large_worker")

# Worker-vs-nymph contrast per species (Stage_Caste groups; several groups
# per side are pooled), used for Hedges' g (section 1) and the cross-species
# comparison (section 2)
caste_contrasts <- list(
  Csec = list(nymph = "Late_Nymph", worker = "Middle_Worker"),
  Mdar = list(nymph = "Late_Nymph", worker = "Middle_Worker"),
  Rgra = list(nymph = "Late_Nymph", worker = "Middle_Worker"),
  Mbar = list(nymph = "Late_Nymph", worker = "Adult_Small_worker")
)

# -----------------------------------------------------------------------------
# 0. One-to-one orthologs and combined cross-species data set
# -----------------------------------------------------------------------------
cat("\n[0] One-to-one orthologs and combined data set\n")
one2one  <- get_one2one_orthologs(ortho_info_demo, species_list = species)
combined <- combine_exp_data(exp_list, one2one, species_list = species)
cat(sprintf("Combined: %d ortholog groups x %d samples\n",
            nrow(combined$abundance), ncol(combined$abundance)))

# -----------------------------------------------------------------------------
# 1. Between-caste divergence: Hedges' g per gene, Worker vs. Nymph (the
#    species-specific groups in `caste_contrasts`). Positive g = higher in
#    workers. Transcriptome-level summary = median |g|.
#
#    The demo abundance data are already corrected for sequencing batch, so
#    the model contains caste group only. For data that are not batch
#    corrected, add batch as a covariate instead, e.g.
#      model = "exp ~ Group + BatchDate"
#    (this needs both groups to be present in at least one shared batch).
# -----------------------------------------------------------------------------
cat("\n[1] Hedges' g (Worker vs. Nymph)\n")

hedges_summary <- list()
hedges_tables  <- list()
for (sp in species) {
  d  <- exp_list[[sp]]
  cc <- caste_contrasts[[sp]]
  keep_samples <- d$sampleInfo$Stage_Caste %in% unlist(cc)
  si <- d$sampleInfo[keep_samples, c("Sample_ID", "Stage_Caste")]
  si$Group <- ifelse(si$Stage_Caste %in% cc$worker, "Worker", "Nymph")
  ab <- d$abundance[, keep_samples]
  # Expressed genes only (median abundance > 1 across the compared samples)
  ab <- ab[apply(ab, 1, median) > 1, ]

  res <- t(sapply(rownames(ab), cal_hedges_g,
                  abundance = ab, sampleInfo = si,
                  model = "exp ~ Group",            # data already batch corrected
                  co_var = "GroupWorker",
                  ref_levels = list(Group = "Nymph")))
  res <- data.frame(gene = rownames(res), res, row.names = NULL)
  hedges_tables[[sp]] <- res

  hedges_summary[[sp]] <- data.frame(
    Species = sp,
    Contrast = paste(paste(cc$worker, collapse = "+"), "vs.", paste(cc$nymph, collapse = "+")),
    n_workers = sum(si$Group == "Worker"),
    n_nymphs = sum(si$Group == "Nymph"),
    n_genes = sum(!is.na(res$Hedges_g)),
    median_abs_g = median(abs(res$Hedges_g), na.rm = TRUE)
  )
}
hedges_summary <- do.call(rbind, hedges_summary)
print(hedges_summary, row.names = FALSE)
write.csv(hedges_summary, file.path(out_dir, "1_hedges_g_summary.csv"),
          row.names = FALSE)

# -----------------------------------------------------------------------------
# 2. Similarity of caste-biased expression between species: angular distance
#    between the species' Worker-vs-Nymph log2 fold-change (LFC) profiles,
#    over one-to-one orthologs that are caste DEGs in at least one species.
#
#    The contrasts differ from section 1 and are set in `lfc_contrasts`.
#
#    NOTE: illustration only. Here, LFC is the difference in mean log2(TPM)
#    between the groups in `lfc_contrasts`, and DEGs are called with limma's
#    moderated t-test (BH-adjusted p < 0.05) on the demo abundance data. In the
#    manuscript, LFCs and DEGs come from formal differential expression
#    analysis; with real data, pass those LFCs to angular_distance_matrix_cols()
#    in the same way.
# -----------------------------------------------------------------------------
cat("\n[2] Angular distance between species' caste-biased expression profiles\n")

# Worker-vs-nymph contrasts for the LFC profiles (Stage_Caste groups)
lfc_contrasts <- list(
  Csec = list(nymph = c("Middle_Nymph", "Late_Nymph"), worker = "Middle_Worker"),
  Mdar = list(nymph = "Early_Nymph",  worker = "Middle_Worker"),
  Rgra = list(nymph = "Early_Nymph",  worker = "Middle_Worker"),
  Mbar = list(nymph = "Early_Nymph",  worker = "Middle_Worker")
)
deg_fdr <- 0.05   # adjusted p-value threshold for calling a caste DEG
lfc_cap <- 3      # |LFC| capped so that a few extreme genes do not dominate

if (!requireNamespace("limma", quietly = TRUE)) stop("Please install 'limma'.")
caste_de <- lapply(setNames(species, species), function(sp) {
  d  <- exp_list[[sp]]
  cc <- lfc_contrasts[[sp]]
  keep <- d$sampleInfo$Stage_Caste %in% c(cc$nymph, cc$worker)
  group <- factor(ifelse(d$sampleInfo$Stage_Caste[keep] %in% cc$worker, "Worker", "Nymph"),
                  levels = c("Nymph", "Worker"))
  ab <- d$abundance[one2one[[sp]], keep, drop = FALSE]
  # Genes constant across the compared samples cannot be tested: LFC = 0, no DEG
  testable <- apply(ab, 1, var) > 0
  out <- data.frame(LFC = rep(0, nrow(ab)), padj = NA_real_, row.names = one2one$GeneID)
  # limma moderated t-test (variance shrinkage across genes, as in count-based
  # DE methods); the Worker coefficient is the LFC in log2(TPM)
  fit <- limma::eBayes(limma::lmFit(ab[testable, ], model.matrix(~ group)))
  tt  <- limma::topTable(fit, coef = "groupWorker", number = Inf, sort.by = "none")
  out$LFC[testable]  <- tt$logFC
  out$padj[testable] <- tt$adj.P.Val
  out
})

lfc  <- sapply(caste_de, `[[`, "LFC");  rownames(lfc)  <- one2one$GeneID
padj <- sapply(caste_de, `[[`, "padj"); rownames(padj) <- one2one$GeneID
is_deg <- !is.na(padj) & padj < deg_fdr
cat("Caste DEGs per species (one-to-one orthologs):\n")
print(data.frame(
  Species  = species,
  Contrast = sapply(lfc_contrasts[species], function(cc)
    paste(cc$worker, "vs.", paste(cc$nymph, collapse = "+"))),
  n_DEG    = colSums(is_deg)[species]
), row.names = FALSE)

# Orthologs that are DEGs in at least one species; LFC capped at +/- lfc_cap
deg_any <- rowSums(is_deg) > 0
lfc_profiles <- pmax(pmin(lfc[deg_any, , drop = FALSE], lfc_cap), -lfc_cap)
cat(sprintf("Orthologs that are DEGs in at least one species: %d\n", nrow(lfc_profiles)))

ang <- angular_distance_matrix_cols(lfc_profiles)
print(round(ang, 3))
write.csv(round(ang, 4), file.path(out_dir, "2_angular_distance_species.csv"))

# Clustering of species by the similarity of their caste-biased expression
pdf(file.path(out_dir, "2_angular_distance_tree.pdf"), width = 5, height = 5)
plot(hclust(as.dist(ang)), hang = -1, main = "Species clustered by caste-bias profiles",
     xlab = "", sub = "", ylab = "Angular distance")
invisible(dev.off())

# Orthologs with the most consistent caste bias: same sign of LFC in all
# species, ranked by the smallest |LFC| among species
same_sign <- abs(rowSums(sign(lfc_profiles))) == length(species)
consistent <- data.frame(
  Orthogroup = rownames(lfc_profiles),
  round(lfc_profiles, 2),
  Bias = ifelse(rowMeans(lfc_profiles) > 0, "Worker", "Nymph"),
  n_species_DEG = rowSums(is_deg[deg_any, , drop = FALSE]),
  min_abs_LFC = round(apply(abs(lfc_profiles), 1, min), 2),
  stringsAsFactors = FALSE
)[same_sign, ]
consistent <- consistent[order(-consistent$min_abs_LFC, -consistent$n_species_DEG), ]

# Annotation: first non-missing value across species, in the order of `species`
annotate_og <- function(og, field) {
  vals <- sapply(species, function(sp) {
    annot <- as.data.frame(gene_info_demo[[sp]])
    gene  <- one2one[[sp]][match(og, one2one$GeneID)]
    annot[[field]][match(gene, annot$GeneID)]
  })
  vals <- matrix(vals, ncol = length(species))
  apply(vals, 1, function(v) v[!is.na(v)][1])
}
consistent$Dmel_name   <- annotate_og(consistent$Orthogroup, "Dmel_name")
consistent$Description <- annotate_og(consistent$Orthogroup, "Description")
write.csv(consistent, file.path(out_dir, "2_consistent_caste_bias.csv"), row.names = FALSE)

cat(sprintf("\nOrthologs biased in the same direction in all species: %d; top 10:\n",
            nrow(consistent)))
top10 <- head(consistent, 10)
top10$Description <- substr(top10$Description, 1, 40)
print(top10[, c("Orthogroup", species, "Bias", "n_species_DEG", "Dmel_name", "Description")],
      row.names = FALSE)

# -----------------------------------------------------------------------------
# 3. Between-caste distance across developmental stages: distance of each
#    nymph (reference) sample to middle-stage workers (target), per species.
# -----------------------------------------------------------------------------
cat("\n[3] Distance of nymph samples to middle-stage workers\n")

stage_order <- c("Early", "Middle", "Late", "Final", "Adult")
dist_all <- list()
for (sp in species) {
  df <- cal_sample_dist(exp_list[[sp]],
                        abundance    = "abundance",
                        target_caste = worker_castes,
                        target_stage = "Middle",
                        ref_caste    = "Nymph",
                        ref_stages   = stage_order,
                        dist_fun     = "man")
  dist_all[[sp]] <- df
}
dist_summary <- do.call(rbind, lapply(dist_all, function(df) {
  df %>% dplyr::group_by(Species, Stage) %>%
    dplyr::summarise(median_distance = median(distance),
              within_worker_baseline = dplyr::first(dist_ref), .groups = "drop")
}))
dist_summary <- as.data.frame(dist_summary)
dist_summary <- dist_summary[order(match(dist_summary$Species, species),
                                   match(dist_summary$Stage, stage_order)), ]
print(dist_summary, row.names = FALSE)
write.csv(dist_summary, file.path(out_dir, "3_sample_distance_summary.csv"),
          row.names = FALSE)

# -----------------------------------------------------------------------------
# 4. Cross-species PCA projection on caste genes: the PC space is built from
#    C. secundus using the top 200 caste DEGs between all nymphs and
#    Middle_Worker (among one-to-one orthologs; Welch t-test on abundance).
#    The other species are projected onto it after ComBat alignment (species
#    as batch, C. secundus kept unchanged). In all species, only early-, middle-
#    and late-stage nymphs (if present) and middle-stage workers are used (`pca_groups`), plus
#    adult workers in Mbar (`pca_groups_extra`).
# -----------------------------------------------------------------------------
cat("\n[4] Reference-based PCA projection on Csec caste DEGs (reference: Csec)\n")

si <- combined$sampleInfo
rownames(si) <- si$Sample_ID

# Caste DEGs in the reference species: all nymphs vs. middle-stage workers
n_deg  <- 200
ref_si <- exp_list[["Csec"]]$sampleInfo
deg_groups <- list(nymph  = unique(ref_si$Stage_Caste[ref_si$Caste == "Nymph"]),
                   worker = "Middle_Worker")
nymph_ids  <- ref_si$Sample_ID[ref_si$Stage_Caste %in% deg_groups$nymph]
worker_ids <- ref_si$Sample_ID[ref_si$Stage_Caste %in% deg_groups$worker]
deg <- t(apply(combined$abundance, 1, function(x) {
  tt <- tryCatch(t.test(x[worker_ids], x[nymph_ids]),
                 error = function(e) NULL)   # e.g. constant expression
  if (is.null(tt)) return(c(diff = NA, t = NA, p = NA))
  c(diff = unname(tt$estimate[1] - tt$estimate[2]),
    t = unname(tt$statistic), p = tt$p.value)
}))
deg <- data.frame(Orthogroup = rownames(deg), deg, row.names = NULL)
deg <- deg[order(deg$p), ]
deg$Bias <- ifelse(deg$diff > 0, "Worker", "Nymph")
top_deg <- head(deg[!is.na(deg$p), ], n_deg)
top_deg$Csec_gene <- one2one$Csec[match(top_deg$Orthogroup, one2one$GeneID)]
annot_csec <- as.data.frame(gene_info_demo[["Csec"]])
top_deg$Dmel_name <- annot_csec$Dmel_name[match(top_deg$Csec_gene, annot_csec$GeneID)]
write.csv(top_deg, file.path(out_dir, "4_csec_caste_DEGs.csv"), row.names = FALSE)
cat(sprintf("Top %d Csec caste DEGs (%s, n = %d vs. %s, n = %d): %d worker-biased, %d nymph-biased; max p = %.2g\n",
            n_deg, deg_groups$worker, length(worker_ids),
            paste(deg_groups$nymph, collapse = "+"), length(nymph_ids),
            sum(top_deg$Bias == "Worker"), sum(top_deg$Bias == "Nymph"),
            max(top_deg$p)))

# Samples used in the PCA (reference and projected species alike)
pca_groups <- c("Early_Nymph", "Middle_Nymph", "Late_Nymph", "Middle_Worker", "Middle_Small_worker")
pca_groups_extra <- list(Mbar = c("Adult_Worker", "Adult_Small_worker"))
groups <- lapply(setNames(species, species), function(sp)
  si$Sample_ID[si$Species == sp &
                 si$Stage_Caste %in% c(pca_groups, pca_groups_extra[[sp]])])
used <- c(pca_groups, unlist(pca_groups_extra))
print(table(si$Species, si$Stage_Caste)[species, intersect(used, unique(si$Stage_Caste))])
pca <- reverse_pca_multi(
  exp_data      = combined$abundance,
  sample.info   = si,
  sample_groups = groups,
  ref_batch     = "Csec",
  genes         = top_deg$Orthogroup,   # PCA on caste DEGs only
  batch_col     = "Species",
  correct_batch = TRUE,
  combat_ref    = TRUE,
  filter        = FALSE,                # genes already selected
  n_pc          = 5
)
cat(sprintf("Projected %d samples onto the PC space of %d Csec samples.\n",
            sum(pca$set_ == "target"), sum(pca$set_ == "ref")))

# Projected samples by species, coloured by caste and shaped by stage
pca_plot <- ggplot2::ggplot(pca, ggplot2::aes(x = PC1, y = PC2, col = Caste, shape = Stage)) +
  ggplot2::geom_point(size = 3) +
  ggplot2::facet_wrap(~Species) +
  ggplot2::theme_bw() +
  ggplot2::labs(title = "Projection onto the Csec caste-DEG PC space")
ggplot2::ggsave(file.path(out_dir, "4_reverse_pca.pdf"), pca_plot, width = 9, height = 7)

# -----------------------------------------------------------------------------
# 5. Illustration: expression of candidate genes across castes and stages in
#    each species: Kr-h1 (Krüppel homolog 1) and Eip93F (Ecdysone-induced
#    protein 93F, E93). One PDF per gene, one page per species.
# -----------------------------------------------------------------------------
candidate_genes <- c("Kr-h1", "Eip93F")
cat("\n[5] Expression plots of", paste(candidate_genes, collapse = " and "), "\n")

for (gene in candidate_genes) {
  pdf(file.path(out_dir, paste0("5_", gene, "_expression.pdf")), width = 7, height = 5)
  for (sp in species) {
    d <- exp_list[[sp]]
    # Order groups by developmental stage for plotting
    sc <- unique(d$sampleInfo[order(match(d$sampleInfo$Stage,
                                          c("Larva", stage_order)),
                                    d$sampleInfo$Caste), "Stage_Caste"])
    d$sampleInfo$Stage_Caste <- factor(d$sampleInfo$Stage_Caste, levels = sc)
    p <- plot_single_candidate_dmel(gene, exp_data = d,
                                    gene_info = as.data.frame(gene_info_demo[[sp]]),
                                    x = "Stage_Caste", col = "Caste",
                                    xlabel = paste(sp, "- stage and caste"))
    print(p)
    cat(sprintf("%s %s: plotted %s\n", gene, sp, p$labels$title))
  }
  invisible(dev.off())
}

cat(sprintf("\nDone. Output written to ./%s/ (run time: %.1f s)\n",
            out_dir, as.numeric(difftime(Sys.time(), t_start, units = "secs"))))
