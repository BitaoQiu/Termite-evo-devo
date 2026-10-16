# =============================================================================
# Functions for comparative analysis of caste-biased expression in termites
#
# 0. Combining species
#      get_one2one_orthologs()        One-to-one orthologs from OrthoFinder output
#      combine_exp_data()             Ortholog-matched multi-species data set
# 1. Between-caste divergence
#      cal_hedges_g()                 Hedges' g per gene (batch-adjusted)
# 2. Similarity of caste-biased expression profiles
#      angular_distance_matrix_cols() Angular distance between profiles
# 3. Between-caste distance across developmental stages
#      cal_dist()                     Sample-to-sample distance wrapper
#      cal_sample_dist()              Target vs. reference sample distances
# 4. Cross-species comparison of caste developmental trajectories
#      reverse_pca_multi()            Reference-based PCA projection
# 5. Illustration
#      plot_single_candidate_dmel()   Expression of one gene across samples
#
# Dependencies: dplyr (CRAN); EDASeq (Bioconductor; used in combine_exp_data);
#               sva (Bioconductor; used in reverse_pca_multi);
#               ggplot2 (CRAN; used in plot_single_candidate_dmel)
# =============================================================================

library(dplyr)


# -----------------------------------------------------------------------------
# 0. Combining species: one-to-one orthologs
# -----------------------------------------------------------------------------

#' One-to-one ortholog table from an OrthoFinder-style orthogroup table
#'
#' Keeps orthogroups with exactly one gene in every species of `species_list`
#' and returns them in the format expected by combine_exp_data().
#'
#' @param ortho_info   Orthogroup table (data frame or list of equal-length
#'                     vectors): a group-ID column plus one column per species,
#'                     holding that species' genes in the group separated by
#'                     `sep` ("" or NA = no gene), as in OrthoFinder's
#'                     Orthogroups.tsv.
#' @param species_list Species that must each have exactly one gene.
#' @param group_col    Name of the group-ID column.
#' @param sep          Separator between gene IDs within a cell (regular expression).
#' @return Data frame with `GeneID` (orthogroup ID) and one gene-ID column per species.
get_one2one_orthologs = function(ortho_info,
                                 species_list,
                                 group_col = "Orthogroup",
                                 sep = ",\\s*") {
  ortho_info <- as.data.frame(ortho_info, stringsAsFactors = FALSE)
  miss <- setdiff(c(group_col, species_list), colnames(ortho_info))
  if (length(miss) > 0) stop("ortho_info lacks column(s): ", paste(miss, collapse = ", "))

  # Number of genes per group and species (empty or NA = 0)
  n_genes <- sapply(species_list, function(sp) {
    x <- trimws(ortho_info[[sp]])
    ifelse(is.na(x) | x == "", 0L, lengths(strsplit(x, sep)))
  })
  n_genes <- matrix(n_genes, ncol = length(species_list))
  keep <- rowSums(n_genes == 1) == length(species_list)

  out <- data.frame(GeneID = ortho_info[[group_col]][keep], stringsAsFactors = FALSE)
  for (sp in species_list) out[[sp]] <- trimws(ortho_info[[sp]][keep])

  # A gene should belong to one group only; groups sharing a gene with another
  # retained group would duplicate its expression, so both are removed
  for (sp in species_list) {
    dup <- out[[sp]] %in% out[[sp]][duplicated(out[[sp]])]
    if (any(dup)) {
      warning(sum(dup), " group(s) removed: ", sp, " gene shared between groups.")
      out <- out[!dup, , drop = FALSE]
    }
  }
  rownames(out) <- NULL
  message(nrow(out), " one-to-one ortholog groups across ",
          paste(species_list, collapse = ", "), ".")
  out
}

#' Combine per-species expression data into one ortholog-matched data set
#'
#' Stacks the expression matrices of several species side by side, with rows
#' matched through a table of one-to-one orthologs, and merges the sample
#' metadata. The output is the input format expected by cross-species
#' analyses such as reverse_pca_multi().
#'
#' Rows are named by the ortholog group ID (`ortho_info$GeneID`). Ortholog
#' groups missing in any species are removed from the abundance matrices
#' (but kept, with NA, in `counts`, if given). The quantile-normalized matrix is computed
#' after this filter.
#'
#' @param exp_data        Named list, one element per species (names = species
#'                        codes, e.g. "Znev"). Each element holds `abundance`
#'                        and `sampleInfo`, and optionally `abundance_correct`
#'                        (used in place of `abundance` when present), `counts`
#'                        and `abundance_weighted_z`. `sampleInfo` needs
#'                        Sample_ID, Caste, Stage and Stage_Caste; Species is
#'                        taken from the list name if absent, and Ind, Colony,
#'                        Tissue, BatchDate and Batch are kept if present
#'                        (NA otherwise).
#' @param ortho_info      One-to-one ortholog table: one row per ortholog group,
#'                        a `GeneID` column with the group ID, and one column per
#'                        species code holding that species' gene ID (as in the
#'                        row names of its expression matrices).
#' @param quantile.method Between-sample normalization passed to
#'                        EDASeq::betweenLaneNormalization(which = ...):
#'                        "upper" (upper-quartile), "median" or "full".
#' @param species_list    Species to combine, in output order. Species absent
#'                        from `exp_data` are dropped with a message.
#' @return List with `abundance`, `abundance.quantile`, `counts` and
#'         `abundance_weighted_z` (each if present in all species), and the merged
#'         `sampleInfo` (one row per column, same order as the matrices).
combine_exp_data = function(exp_data,
                            ortho_info,
                            quantile.method = "upper",
                            species_list = c("Mdar", "Mvia", "Znev", "Csec",
                                             "Pcan", "Rgra", "Mbar")) {
  if (!requireNamespace("EDASeq", quietly = TRUE)) stop("Please install 'EDASeq'.")

  # --- Check inputs -----------------------------------------------------------
  missing_sp <- setdiff(species_list, names(exp_data))
  if (length(missing_sp) > 0) {
    message("Not in exp_data, skipped: ", paste(missing_sp, collapse = ", "))
  }
  species_list <- intersect(species_list, names(exp_data))
  if (length(species_list) == 0) stop("None of species_list is in exp_data.")

  if (!"GeneID" %in% colnames(ortho_info)) stop("ortho_info needs a 'GeneID' column.")
  no_ortho_col <- setdiff(species_list, colnames(ortho_info))
  if (length(no_ortho_col) > 0) {
    stop("ortho_info has no column for: ", paste(no_ortho_col, collapse = ", "))
  }

  # Matrix columns and metadata rows must refer to the same samples in the same
  # order, because the merged sampleInfo is stacked in the same way as the matrices
  for (sp in species_list) {
    if (!identical(colnames(exp_data[[sp]]$abundance),
                   as.character(exp_data[[sp]]$sampleInfo$Sample_ID))) {
      stop("In ", sp, ", columns of 'abundance' do not match sampleInfo$Sample_ID.")
    }
  }

  # --- Stack expression matrices by ortholog group ----------------------------
  # For each species, pick the rows of its orthologs in the order of ortho_info
  # (NA rows where a species lacks the gene), then bind species column-wise
  stack_orthologs <- function(slot) {
    mat <- do.call(cbind, lapply(species_list, function(sp) {
      m <- exp_data[[sp]][[slot]]
      matched_idx <- match(ortho_info[[sp]], rownames(m))
      m[matched_idx, , drop = FALSE]
    }))
    rownames(mat) <- ortho_info$GeneID
    mat
  }
  drop_incomplete <- function(mat) mat[!apply(mat, 1, anyNA), , drop = FALSE]

  # Batch-corrected abundance replaces the raw abundance where available
  for (sp in species_list) {
    if (!is.null(exp_data[[sp]]$abundance_correct)) {
      exp_data[[sp]]$abundance <- exp_data[[sp]]$abundance_correct
    }
  }

  # Optional matrices: combined only if every species has them
  in_all_species <- function(slot) {
    all(vapply(species_list, function(sp) !is.null(exp_data[[sp]][[slot]]), logical(1)))
  }

  exp_data_full <- list()
  if (in_all_species("counts")) {
    exp_data_full$counts <- stack_orthologs("counts")
  }
  if (in_all_species("abundance_weighted_z")) {
    exp_data_full$abundance_weighted_z <- drop_incomplete(stack_orthologs("abundance_weighted_z"))
  }

  exp_data_full$abundance <- drop_incomplete(stack_orthologs("abundance"))
  if (nrow(exp_data_full$abundance) == 0) {
    stop("No ortholog group is present in all species.")
  }

  # Between-sample normalization across all species (e.g. upper-quartile), so
  # that library-wide differences between species do not dominate comparisons
  exp_data_full$abundance.quantile <- EDASeq::betweenLaneNormalization(
    exp_data_full$abundance, which = quantile.method, offset = FALSE, round = FALSE)

  # --- Harmonize species-specific metadata ------------------------------------
  # Bger (cockroach, non-social outgroup): stage labels shortened, and all
  # individuals treated as the "Nymph" developmental path
  if ("Bger" %in% species_list) {
    si <- exp_data[["Bger"]]$sampleInfo
    si$Stage <- as.character(factor(
      si$Stage,
      levels = c("Egg day 13", "Larva", "Early", "Middle", "Late", "Adult"),
      labels = c("Egg", "Larva", "Early", "Middle", "Late", "Adult")))
    si$Caste <- "Nymph"
    si$Stage_Caste <- paste(si$Stage, "Nymph", sep = "_")
    exp_data[["Bger"]]$sampleInfo <- si
  }
  # Zang: worker samples only, stages up to "Middle" (other labels become NA);
  # no colony information
  if ("Zang" %in% species_list) {
    si <- exp_data[["Zang"]]$sampleInfo
    si$Stage <- as.character(factor(si$Stage,
                                    levels = c("Larva", "Early", "Early-Middle", "Middle")))
    si$Stage_Caste <- paste(si$Stage, "Worker", sep = "_")
    si$Colony <- NA
    exp_data[["Zang"]]$sampleInfo <- si
  }

  # --- Merge sample metadata --------------------------------------------------
  target_atr <- c("Sample_ID", "Ind", "Caste", "Colony", "Stage", "Stage_Caste",
                  "Tissue", "Species", "BatchDate", "Batch")
  required_atr <- c("Sample_ID", "Caste", "Stage", "Stage_Caste")
  for (sp in species_list) {
    si <- exp_data[[sp]]$sampleInfo
    miss <- setdiff(required_atr, colnames(si))
    if (length(miss) > 0) {
      stop("sampleInfo of ", sp, " lacks column(s): ", paste(miss, collapse = ", "))
    }
    # Species from the list name if absent; other bookkeeping columns NA
    if (is.null(si$Species)) si$Species <- sp
    for (col in setdiff(target_atr, colnames(si))) si[[col]] <- NA
    exp_data[[sp]]$sampleInfo <- si
  }
  exp_data_full$sampleInfo <- do.call(rbind, lapply(species_list, function(sp) {
    exp_data[[sp]]$sampleInfo[, target_atr]
  }))

  # Factor levels fix the plotting/ordering of stages and castes.
  # Labels not in these lists would become NA, so they are reported.
  stage_levels <- c("Egg", "Larva", "Early", "Early-Middle", "Middle", "Late",
                    "Final", "Adult", "Adult_Late")
  caste_levels <- c("Queen", "King", "Queen/King", "Major worker", "Nymph",
                    "Worker", "Large_worker", "Small_worker",
                    "Small_soldier", "Large_soldier")
  check_levels <- function(x, levels, what) {
    unknown <- setdiff(unique(na.omit(x)), levels)
    if (length(unknown) > 0) {
      warning(what, " label(s) not in the level list, set to NA: ",
              paste(unknown, collapse = ", "))
    }
  }
  check_levels(exp_data_full$sampleInfo$Stage, stage_levels, "Stage")
  check_levels(exp_data_full$sampleInfo$Caste, caste_levels, "Caste")

  si <- exp_data_full$sampleInfo
  si$Stage <- factor(si$Stage, levels = stage_levels)

  # Stage_raw: stage scheme harmonized across species. In Csec, Pcan and Znev,
  # "Late" samples are grouped with "Middle"
  si$Stage_raw <- si$Stage
  si$Stage_raw[si$Stage == "Late" & si$Species %in% c("Csec", "Pcan", "Znev")] <- "Middle"

  si$Caste   <- factor(si$Caste, levels = caste_levels)
  si$Species <- factor(si$Species, levels = species_list)
  si$Species_Stage_Caste <- paste(si$Species, si$Stage, si$Caste, sep = "_")

  # Worker type of each species (FW / TW / SW); species not listed are FW
  si$Worker_Type <- "FW"
  si$Worker_Type[si$Species %in% c("Mdar", "Rgra")] <- "TW"
  si$Worker_Type[si$Species %in% c("Mvia", "Mbar", "Odon", "Trin", "Mbel")] <- "SW"
  si$Worker_Type <- factor(si$Worker_Type, levels = c("FW", "TW", "SW"))

  exp_data_full$sampleInfo <- si
  exp_data_full
}


# -----------------------------------------------------------------------------
# 1. Between-caste divergence: Hedges' g
# -----------------------------------------------------------------------------

#' Hedges' g for one gene from a batch-adjusted linear model
#'
#' Quantifies between-caste divergence for a single gene as a standardized
#' effect size. For a transcriptome-level summary, take the median across genes.
#'
#' @param geneID     Row name of the gene in `abundance`.
#' @param abundance  Expression matrix (genes x samples).
#' @param sampleInfo Sample metadata (one row per sample).
#' @param model      Model formula as a string; the response must be named `exp`.
#' @param co_var     Coefficient to convert, as named in coef(fit),
#'                   e.g. "CasteWorker" = Worker vs. reference caste.
#' @param ref_levels Reference level per factor; effects are relative to these.
#' @return coef, residual SD, Cohen's d, Hedges' g (NA if not estimable).
cal_hedges_g = function(geneID,
                        abundance,
                        sampleInfo,
                        model = "exp ~ Caste + BatchDate",  # adjusts for batch
                        co_var = "CasteWorker",
                        ref_levels = list(Caste = "Nymph")) {
  # Assumes columns of `abundance` match rows of `sampleInfo` in order
  tmp_data <- cbind(
    exp = as.numeric(abundance[geneID, ]),
    sampleInfo,
    stringsAsFactors = FALSE
  )
  tmp_data <- droplevels(tmp_data)
  for (var in names(ref_levels)) {
    tmp_data[[var]] <- relevel(factor(tmp_data[[var]]), ref = ref_levels[[var]])
  }
  fit <- lm(as.formula(model), data = tmp_data)
  
  coef_est <- unname(coef(fit)[co_var])
  # Standardizer: residual SD after removing caste and batch effects
  # (used in place of the raw pooled SD of classic Cohen's d)
  rsd <- sigma(fit)
  
  # NA coef = contrast not estimable (e.g. caste confounded with batch)
  if (is.na(coef_est) || is.na(rsd) || rsd == 0) {
    return(c(coef = NA_real_, Rsd = rsd, Cohens_d = NA_real_, Hedges_g = NA_real_))
  }
  
  d <- coef_est / rsd
  # Hedges' small-sample correction, with residual df in place of n1 + n2 - 2
  df_res <- df.residual(fit)
  J <- 1 - 3 / (4 * df_res - 1)
  g <- d * J
  c(coef = coef_est, Rsd = rsd, Cohens_d = d, Hedges_g = g)
}


# -----------------------------------------------------------------------------
# 2. Similarity of caste-biased expression profiles: angular distance
# -----------------------------------------------------------------------------

#' Angular distance between all pairs of columns
#' (e.g. log fold-change profiles of caste expression differences)
#'
#' Angular distance = acos(cosine similarity) / pi, scaled to [0, 1]:
#' 0 = same direction, 0.5 = orthogonal, 1 = opposite.
#' Unlike 1 - cosine, it is a proper metric (satisfies the triangle inequality).
#' Only the direction of each profile matters, not its magnitude.
#' For non-negative input (e.g. raw abundance) the range is [0, 0.5].
#'
#' @param mat Numeric matrix; distances are computed between columns.
#' @return Symmetric columns x columns distance matrix.
angular_distance_matrix_cols = function(mat) {
  X <- as.matrix(mat)
  # Missing values are set to 0 (= no change, for fold-change input)
  X[is.na(X)] <- 0
  storage.mode(X) <- "double"
  
  # Cosine similarity for all column pairs at once: X'X / (||x_i|| * ||x_j||)
  G <- crossprod(X)
  nrm <- sqrt(colSums(X * X))
  den <- outer(nrm, nrm, "*")
  
  cosine <- as.numeric(G) / as.numeric(den)
  dim(cosine) <- dim(G)
  
  # All-zero columns give 0/0 (NaN). They are treated as orthogonal to all other
  # columns (distance 0.5) but identical to themselves (distance 0).
  # This must happen before clamping and acos().
  cosine[!is.finite(cosine)] <- 0
  zero_cols <- which(nrm == 0)
  if (length(zero_cols) > 0) {
    diag(cosine)[zero_cols] <- 1
  }
  
  # Floating-point error can push values slightly outside [-1, 1],
  # which would make acos() return NaN
  cosine <- pmax(-1, pmin(1, cosine))
  dist_vec <- acos(cosine) / pi
  # pmin/pmax drop the matrix dimensions, so the matrix is rebuilt here
  dist_mat <- matrix(dist_vec, nrow = nrow(G), ncol = ncol(G))
  
  cn <- colnames(X)
  if (!is.null(cn)) dimnames(dist_mat) <- list(cn, cn)
  dist_mat
}


# -----------------------------------------------------------------------------
# 3. Between-caste transcriptomic distance across developmental stages
# -----------------------------------------------------------------------------

#' Pairwise distances between columns (samples)
#'
#' "s"/"p"/"k": 1 - Spearman/Pearson/Kendall correlation (complete observations);
#' "angular": angular distance; anything else is passed to dist(), e.g. "man".
#'
#' @param tmp.abundance Expression matrix (genes x samples).
#' @param dist_fun      Distance code (see above).
#' @return A `dist` object.
cal_dist = function(tmp.abundance, dist_fun) {
  if (dist_fun %in% c("s", "p", "k")) {
    tmp.Dists = as.dist(1 - cor(tmp.abundance, method = dist_fun, use = "c"))
  } else if (dist_fun == "angular") {
    tmp.Dists = as.dist(angular_distance_matrix_cols(tmp.abundance))
  } else {
    tmp.Dists = dist(t(tmp.abundance), method = dist_fun)
  }
  return(tmp.Dists)
}

#' Expression distances between target-caste samples and reference samples
#'
#' Compares target samples (e.g. mid-stage workers) to reference samples
#' (e.g. nymphs at each stage) to show which reference stage the targets
#' most resemble.
#'
#' @param exp_data         List with `sampleInfo` and expression matrices.
#' @param abundance        Name of the expression matrix inside `exp_data`.
#' @param subset_BatchDate Optional batch(es) to restrict the analysis to.
#' @param dist_fun         Distance code passed to cal_dist().
#' @return Long data frame of target x reference distances, with sample metadata,
#'         the per-reference median distance, and a within-target baseline (dist_ref).
cal_sample_dist = function(
    exp_data,
    abundance = "abundance",
    target_caste = c("Small_worker", "Worker"),
    target_stage = "Middle",
    ref_caste = "Nymph",
    ref_stages = c("Early", "Middle", "Late", "Adult"),
    subset_BatchDate = NULL,
    dist_fun = "man") {
  if (!is.null(subset_BatchDate)) {
    subset_samples <- exp_data$sampleInfo$BatchDate %in% subset_BatchDate
  } else {
    subset_samples <- rep(TRUE, nrow(exp_data$sampleInfo))
  }
  sample_info <- exp_data$sampleInfo[subset_samples, , drop = FALSE]
  abundance_mat <- exp_data[[abundance]][, subset_samples, drop = FALSE]
  
  # Keep expressed genes only (median abundance > 1 across retained samples),
  # so that noise from near-zero genes does not dominate the distances
  keep_rows <- apply(abundance_mat, 1, median) > 1
  abundance_filtered <- abundance_mat[keep_rows, , drop = FALSE]
  
  target_samples <- sample_info$Sample_ID[
    sample_info$Stage %in% target_stage &
      sample_info$Caste %in% target_caste
  ]
  
  ref_samples <- sample_info$Sample_ID[
    sample_info$Stage %in% ref_stages &
      sample_info$Caste %in% ref_caste
  ]
  
  # Assumes columns of the abundance matrix are named by Sample_ID
  dist_full <- as.matrix(cal_dist(abundance_filtered, dist_fun = dist_fun))
  # Correlation-based distances (s/p/k) and angular distance are bounded and
  # left as is. Other distances (e.g. Manhattan) grow with the number of genes,
  # so they are scaled to a per-gene average to stay comparable across
  # gene filters and batch subsets.
  if (!dist_fun %in% c("s", "p", "k", "angular")) {
    dist_full <- dist_full / nrow(abundance_filtered)
  }
  dist_subset <- dist_full[target_samples, ref_samples, drop = FALSE]
  
  # Baseline: average pairwise distance among target samples themselves,
  # i.e. how close "close" is. Upper triangle only, which excludes self-distances
  # (diagonal zeros) and counts each pair once.
  tt <- dist_full[target_samples, target_samples, drop = FALSE]
  dist_ref <- mean(tt[upper.tri(tt)], na.rm = TRUE)
  
  df_long <- as.data.frame(as.table(dist_subset), stringsAsFactors = FALSE)
  colnames(df_long) <- c("row_sample", "col_sample", "distance")
  
  # Rows = target samples, columns = reference samples
  row_info <- sample_info[, c("Sample_ID", "Caste", "Stage")]
  colnames(row_info) <- c("row_sample", "row_Caste", "row_Stage")
  
  col_info <- sample_info[, c("Sample_ID", "Caste", "Stage", "Species")]
  colnames(col_info) <- c("col_sample", "Caste", "Stage", "Species")
  
  df_long2 <- df_long %>%
    dplyr::left_join(row_info, by = "row_sample") %>%
    dplyr::left_join(col_info, by = "col_sample") %>%
    # Summary distance of each reference sample to all target samples
    # (median, robust to single outlier targets)
    dplyr::group_by(col_sample) %>%
    dplyr::mutate(median_dist = median(distance, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    # Order reference samples by stage, then by closeness to targets (for plotting)
    dplyr::arrange(Stage, median_dist) %>%
    dplyr::mutate(
      col_sample = factor(col_sample, levels = unique(col_sample)),
      dist_ref = dist_ref
    )
  
  df_long2
}


# -----------------------------------------------------------------------------
# 4. Principal component projection: comparing caste developmental
#    trajectories across species
# -----------------------------------------------------------------------------

#' Reference-based ("reverse") PCA across species
#'
#' Builds the PC space from the reference group only (e.g. one species), then
#' projects the target groups (other species) into that space. The PCs therefore
#' describe the major axes of variation in the reference, and the targets are
#' placed along those same axes without influencing them.
#' Centering and scaling also use reference statistics only.
#'
#' @param exp_data      Expression matrix (genes x samples).
#' @param sample.info   Sample metadata (rows = samples).
#' @param sample_groups Named list of sample IDs or column indices per group.
#' @param ref_batch     Group that defines the PC space.
#' @param correct_batch If TRUE, apply ComBat across `batch_col` before PCA.
#' @param combat_ref    If TRUE, ComBat uses the reference as `ref.batch`: reference
#'                      data stay unchanged and targets are mapped onto it.
#' @param filter        If TRUE, keep the `n_gene` most variable genes in the reference.
#' @return Data frame of PC scores (PC1..PCn) for reference and target samples,
#'         joined with sample metadata; `set_` marks "ref" vs "target".
reverse_pca_multi <- function(
    exp_data,
    sample.info,
    sample_groups,            # named list: list(Csec = idx_or_ids, Pcan = idx_or_ids, ...)
    ref_batch,                # e.g. "Csec"
    target_batches = NULL,    # optional character vector, e.g. c("Pcan","Znev")
    genes = NULL,
    sample_id_col = NULL,     # set if rownames(sample.info) are not sample IDs
    batch_col = "Species",    # generic metadata column used for ComBat if available
    correct_batch = TRUE,
    combat_ref = TRUE,        # if TRUE, ComBat keeps the reference unchanged (ref.batch)
    center = TRUE,
    scale = FALSE,
    method_d = "combat",
    filter = FALSE,
    n_gene = 512,
    n_pc = 5
) {
  if (is.null(names(sample_groups)) || any(names(sample_groups) == "")) {
    stop("sample_groups must be a named list.")
  }
  if (!ref_batch %in% names(sample_groups)) {
    stop("ref_batch is not a name in sample_groups.")
  }
  
  # Default: all groups other than the reference are projected
  if (is.null(target_batches)) {
    target_batches <- setdiff(names(sample_groups), ref_batch)
  } else {
    if (!all(target_batches %in% names(sample_groups))) {
      stop("Some target_batches are not names in sample_groups.")
    }
    target_batches <- setdiff(target_batches, ref_batch)
  }
  if (length(target_batches) == 0) {
    stop("No valid target_batches after excluding ref_batch.")
  }
  
  if (!is.null(sample_id_col)) {
    rownames(sample.info) <- as.character(sample.info[[sample_id_col]])
  }
  
  if (is.null(genes)) genes <- rownames(exp_data)
  genes <- intersect(genes, rownames(exp_data))
  if (length(genes) == 0) stop("No valid genes found in exp_data.")
  
  # Align metadata rows to expression columns by sample ID
  exp_data <- exp_data[genes, , drop = FALSE]
  all_samples <- colnames(exp_data)
  sample.info <- sample.info[all_samples, , drop = FALSE]
  
  # Groups may be given as column indices or as sample IDs
  to_ids <- function(x) {
    if (is.numeric(x)) all_samples[x] else as.character(x)
  }
  
  used_groups <- c(ref_batch, target_batches)
  group_ids <- lapply(sample_groups[used_groups], to_ids)
  
  all_ids <- unlist(group_ids, use.names = FALSE)
  dup_ids <- unique(all_ids[duplicated(all_ids)])
  if (length(dup_ids) > 0) {
    stop("Some samples are assigned to multiple groups: ", paste(dup_ids, collapse = ", "))
  }
  
  # Reference samples define the PC space ("train");
  # target samples are only projected ("test")
  train_ids <- unique(group_ids[[ref_batch]])
  test_ids <- unique(unlist(group_ids[target_batches], use.names = FALSE))
  
  if (length(train_ids) == 0) stop("No training samples in ref_batch.")
  if (length(test_ids) == 0) stop("No testing samples in target_batches.")
  
  use_ids <- unique(c(train_ids, test_ids))
  exp_use <- exp_data[, use_ids, drop = FALSE]
  # Missing values are treated as zero expression
  exp_use[is.na(exp_use)] <- 0
  
  info_use <- sample.info[use_ids, , drop = FALSE]
  
  sample_to_group <- setNames(rep(NA_character_, length(use_ids)), use_ids)
  for (g in names(group_ids)) {
    sample_to_group[group_ids[[g]]] <- g
  }
  info_use$list_batch_ <- sample_to_group[rownames(info_use)]
  info_use$set_ <- ifelse(info_use$list_batch_ == ref_batch, "ref", "target")
  
  # Optional gene filter: top-variance genes, ranked in the reference only,
  # so that target samples do not influence gene selection
  if (filter) {
    ref_var <- apply(exp_use[, train_ids, drop = FALSE], 1, var, na.rm = TRUE)
    n_use <- min(n_gene, length(ref_var))
    keep_genes <- names(sort(ref_var, decreasing = TRUE))[seq_len(n_use)]
    exp_use <- exp_use[keep_genes, , drop = FALSE]
  }
  
  # ComBat removes per-gene location/scale differences between batches
  # (by default, species). mod = ~1: no biological covariate is protected,
  # so any group difference aligned with batch is removed as well.
  # With combat_ref = TRUE, the reference batch is left unchanged and targets
  # are mapped onto it, so the PC space depends on the reference only.
  # With FALSE, all batches (including the reference) are mapped to a pooled
  # average, which lets target samples influence the reference data.
  if (correct_batch) {
    if (method_d != "combat") stop("Only method_d='combat' is implemented.")
    if (!requireNamespace("sva", quietly = TRUE)) stop("Please install 'sva'.")
    
    # Falls back to the sample_groups labels if batch_col is not in the metadata
    if (!is.null(batch_col) && batch_col %in% colnames(info_use)) {
      batch <- factor(info_use[[batch_col]])
    } else {
      batch <- factor(info_use$list_batch_)
    }
    
    ref_level <- NULL
    if (combat_ref) {
      # ComBat needs the reference as a level of `batch`, which may be labelled
      # differently from ref_batch (e.g. species name vs. group name),
      # so it is taken from the reference samples themselves
      ref_level <- unique(as.character(batch[info_use$set_ == "ref"]))
      if (length(ref_level) != 1) {
        stop("Reference samples span multiple batch levels in '", batch_col,
             "': ", paste(ref_level, collapse = ", "),
             ". ref.batch requires a single reference batch.")
      }
      # Targets sharing the reference batch label would be left unadjusted
      if (any(batch[info_use$set_ == "target"] == ref_level)) {
        warning("Some target samples share the reference batch label '",
                ref_level, "' and will not be adjusted by ComBat.")
      }
    }
    
    modcombat <- model.matrix(~1, data = info_use)
    exp_use <- sva::ComBat(
      dat = exp_use, batch = batch, mod = modcombat,
      mean.only = FALSE, par.prior = TRUE, prior.plots = FALSE,
      ref.batch = ref_level
    )
  }
  
  train_mat <- exp_use[, train_ids, drop = FALSE]
  test_mat  <- exp_use[, test_ids, drop = FALSE]
  
  # Center (and optionally scale) both sets with reference statistics,
  # so targets are expressed relative to the reference
  if (center) {
    ref_mean <- rowMeans(train_mat)
    train_center <- sweep(train_mat, 1, ref_mean, "-")
    test_center  <- sweep(test_mat, 1, ref_mean, "-")
  } else {
    train_center <- train_mat
    test_center  <- test_mat
  }
  if (scale) {
    ref_sd <- apply(train_center, 1, sd)
    # Genes constant in the reference: leave unscaled to avoid division by zero
    ref_sd[ref_sd == 0 | is.na(ref_sd)] <- 1
    train_center <- sweep(train_center, 1, ref_sd, "/")
    test_center  <- sweep(test_center, 1, ref_sd, "/")
  }
  n_pc <- min(n_pc, ncol(train_center), nrow(train_center))
  # PCA via SVD of the reference data; columns of u are the gene loadings
  sv <- svd(train_center, nu = n_pc, nv = n_pc)
  
  # PC scores = projection onto reference loadings (same formula for both sets).
  # PC signs are arbitrary, as in any PCA.
  pc_train <- t(t(sv$u) %*% train_center)
  pc_test  <- t(t(sv$u) %*% test_center)
  
  out <- rbind(pc_train, pc_test)
  out <- data.frame(out, info_use[rownames(out), , drop = FALSE], check.names = FALSE)
  colnames(out)[seq_len(n_pc)] <- paste0("PC", seq_len(n_pc))
  
  return(out)
}

# -----------------------------------------------------------------------------
# 5. Illustration: expression of a single candidate gene
# -----------------------------------------------------------------------------

#' Plot the expression of one gene, selected by its Drosophila ortholog name
#'
#' Several termite genes can share the same Drosophila best-hit name; the one
#' with the highest bit score is plotted.
#'
#' @param dmel_name  Drosophila gene name (as in gene_info$Dmel_name), e.g. "Kr-h1".
#' @param exp_data   One species' data: list with an expression matrix (genes x
#'                   samples) and `sampleInfo` in the same sample order.
#' @param gene_info  Gene annotation of that species, with columns GeneID,
#'                   Dmel_name, Dmel_bitScore and Description.
#' @param x, col     sampleInfo columns for the x axis and the colour.
#' @param shape      sampleInfo column for the point shape, or a fixed shape code.
#' @param xlabel     x-axis label.
#' @param plot_type  "boxplot" (box plot with samples overlaid) or "point".
#' @param abundance  Name of the expression matrix inside `exp_data`.
#' @param size       Point size.
#' @return A ggplot object.
plot_single_candidate_dmel = function(dmel_name,
                                      exp_data,
                                      gene_info,
                                      shape = 20,
                                      x = "Caste",
                                      col = "Caste",
                                      xlabel = "Phenotype",
                                      plot_type = c("boxplot", "point"),
                                      abundance = "abundance",
                                      size = 2) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Please install 'ggplot2'.")
  plot_type <- match.arg(plot_type)
  
  # Genes with this Drosophila name that are present in the expression matrix;
  # the best hit (highest bit score) is used
  mat <- exp_data[[abundance]]
  gene_candidates <- gene_info[gene_info$Dmel_name %in% dmel_name &
                                 gene_info$GeneID %in% rownames(mat), , drop = FALSE]
  if (nrow(gene_candidates) == 0) {
    stop("No gene with Dmel_name '", dmel_name, "' found in the expression data.")
  }
  best <- order(gene_candidates$Dmel_bitScore, decreasing = TRUE)[1]
  geneid    <- gene_candidates$GeneID[best]
  gene_desc <- gene_candidates$Description[best]
  
  gene.data <- data.frame(Expression.level = as.numeric(mat[geneid, ]),
                          exp_data$sampleInfo,
                          check.names = FALSE)
  
  # Shape: a sampleInfo column if one is named, otherwise a fixed point shape
  shape_aes <- if (is.character(shape) && shape %in% colnames(gene.data)) {
    ggplot2::aes(shape = .data[[shape]])
  } else {
    ggplot2::aes(shape = I(shape))
  }
  
  p1 <- ggplot2::ggplot(gene.data,
                        ggplot2::aes(x = .data[[x]], colour = .data[[col]],
                                     y = .data[["Expression.level"]])) +
    shape_aes +
    ggplot2::theme_bw() +
    ggplot2::ylab("Expression abundance: log2(TPM)") +
    ggplot2::xlab(xlabel) +
    ggplot2::labs(title = geneid,
                  subtitle = if (is.na(gene_desc) || gene_desc == "") dmel_name
                             else paste0(dmel_name, " (", gene_desc, ")")) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  
  if (plot_type == "boxplot") {
    # Outliers are not drawn by the box plot, since every sample is shown as a point
    p1 <- p1 +
      ggplot2::geom_boxplot(position = ggplot2::position_dodge(0.8), outlier.shape = NA) +
      ggplot2::geom_point(position = ggplot2::position_jitterdodge(
        jitter.width = 0.2, dodge.width = 0.8, seed = 1),
        size = size, alpha = 0.5)
  } else {
    p1 <- p1 + ggplot2::geom_point(size = size)
  }
  p1
}
