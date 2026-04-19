#!/usr/bin/env Rscript
# =============================================================================
# calibrate_mcra_thresholds.R
#
# Empirically calibrates phylogenetic distance thresholds for mcrA at each
# taxonomic rank (species → phylum) using the Youden-J (ROC) method.
#
# Method (mirrors Yarza et al. 2014, Nat Rev Microbiol 12:635):
#   1. Load sequences and taxonomy from mcrAtemplate.fasta + tax4mcrA.taxonomy
#   2. Align with MAFFT; build GTR tree with FastTree (NJ fallback)
#   3. Compute full patristic distance matrix
#   4. For every sequence pair, record the *deepest rank* at which they differ
#   5. For each rank R, treat it as a binary classification problem:
#        positive  = pairs that differ AT rank R (inter-rank boundary)
#        negative  = pairs that differ BELOW rank R (within-rank)
#      Find the K80 distance cutoff d* that maximises the Youden index
#        J(d) = sensitivity(d) + specificity(d) - 1
#   6. Report d* per rank with confidence intervals (bootstrap, n = 999)
#
# Output:
#   mcra_rank_thresholds.tsv  — thresholds table (copy into build_mcrA_db.R)
#   mcra_rank_thresholds.pdf  — ROC curves and distance distributions per rank
#
# Usage:
#   Rscript calibrate_mcra_thresholds.R
#   (run from the folder containing mcrAtemplate.fasta and tax4mcrA.taxonomy)
# =============================================================================

# ── Dependencies ─────────────────────────────────────────────────────────────
required_pkgs <- c("ape", "phangorn", "ggplot2", "dplyr", "tidyr", "pROC", "gridExtra")
missing_pkgs  <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0L)
  stop("Install missing packages first:\n  install.packages(c(",
       paste(sprintf('"%s"', missing_pkgs), collapse = ", "), "))")

suppressPackageStartupMessages({
  library(ape)
  library(phangorn)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(pROC)
  library(gridExtra)
})

set.seed(42)

# ── Configuration ─────────────────────────────────────────────────────────────
FASTA_FILE      <- "mcrAtemplate.fasta"
TAX_FILE        <- "tax4mcrA.taxonomy"
FASTTREE_PATH   <- "~/anaconda3/bin/FastTree"
MAFFT_PATH      <- "mafft"
MAFFT_THREADS   <- 8L
N_BOOT          <- 999L          # bootstrap replicates for CI
MIN_PAIRS       <- 30L           # minimum pairs per rank for reliable threshold
OUT_TSV         <- "mcra_rank_thresholds.tsv"
OUT_PDF         <- "mcra_rank_thresholds.pdf"

TAX_LEVELS <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# Ranks we calibrate (species needs many intra-species pairs; often sparse)
CALIBRATE_RANKS <- c("Species", "Genus", "Family", "Order", "Class", "Phylum")

# ── Helper: expand ~ paths ────────────────────────────────────────────────────
resolve_bin <- function(name) {
  exp <- path.expand(name)
  if (file.exists(exp) && file.access(exp, 1L) == 0L) return(exp)
  found <- Sys.which(name)
  if (nzchar(found)) return(found)
  lc <- tolower(name)
  found2 <- Sys.which(lc)
  if (nzchar(found2)) return(found2)
  name
}

# ── Helper: read FASTA ────────────────────────────────────────────────────────
read_fasta <- function(path) {
  lines <- readLines(path)
  hdr_idx <- grep("^>", lines)
  starts  <- hdr_idx + 1L
  ends    <- c(hdr_idx[-1L] - 1L, length(lines))
  seqs    <- mapply(function(s, e) paste(lines[s:e], collapse = ""),
                    starts, ends, SIMPLIFY = TRUE)
  names(seqs) <- sub("^>", "", lines[hdr_idx])
  seqs
}

# =============================================================================
# 1. LOAD DATA
# =============================================================================
message("\n── 1. Loading sequences and taxonomy ──────────────────────────────────")

if (!file.exists(FASTA_FILE)) stop("FASTA not found: ", FASTA_FILE)
if (!file.exists(TAX_FILE))   stop("Taxonomy not found: ", TAX_FILE)

raw_seqs <- read_fasta(FASTA_FILE)
message("  Sequences loaded: ", length(raw_seqs))

# Parse taxonomy file → named list: accession → named character vector
tax_lines <- readLines(TAX_FILE)
tax_map   <- list()
for (ln in tax_lines) {
  ln <- trimws(ln)
  if (!nzchar(ln)) next
  parts <- strsplit(ln, "\t")[[1L]]
  if (length(parts) < 2L) next
  acc    <- trimws(parts[1L])
  levels <- trimws(strsplit(gsub(";$", "", trimws(parts[2L])), ";")[[1L]])
  while (length(levels) < length(TAX_LEVELS))
    levels <- c(levels, "Unclassified")
  tax_map[[acc]] <- setNames(levels[seq_along(TAX_LEVELS)], TAX_LEVELS)
}
message("  Taxonomy entries : ", length(tax_map))

# Keep sequences that have taxonomy AND are not fully Unclassified
keep   <- intersect(names(raw_seqs), names(tax_map))
seqs   <- raw_seqs[keep]
tax_df <- do.call(rbind, lapply(keep, function(a) {
  v <- tax_map[[a]]
  as.data.frame(t(v), stringsAsFactors = FALSE)
}))
rownames(tax_df) <- keep

# Filter: require at least Kingdom to be classified
classified <- rownames(tax_df)[tax_df$Kingdom != "Unclassified"]
seqs   <- seqs[classified]
tax_df <- tax_df[classified, ]
message("  Sequences with taxonomy: ", length(seqs))

if (length(seqs) < 10L)
  stop("Too few classified sequences to calibrate thresholds.")

# =============================================================================
# 2. ALIGN
# =============================================================================
message("\n── 2. Aligning with MAFFT ──────────────────────────────────────────────")

tmp_dir    <- tempdir()
input_fa   <- file.path(tmp_dir, "calib_input.fasta")
aligned_fa <- file.path(tmp_dir, "calib_aligned.fasta")

# Write input
writeLines(
  unlist(mapply(function(nm, sq) c(paste0(">", nm), sq),
                names(seqs), seqs, SIMPLIFY = FALSE)),
  input_fa
)

mafft_bin <- resolve_bin(MAFFT_PATH)
cmd_aln <- sprintf("%s --auto --thread %d --quiet %s > %s 2>/dev/null",
                   shQuote(mafft_bin), MAFFT_THREADS,
                   shQuote(input_fa), shQuote(aligned_fa))
aln_ret <- system(cmd_aln)
if (aln_ret != 0L || !file.exists(aligned_fa) || file.info(aligned_fa)$size == 0L)
  stop("MAFFT alignment failed. Check that mafft is installed.")
message("  Alignment done: ", aligned_fa)

# =============================================================================
# 3. BUILD TREE
# =============================================================================
message("\n── 3. Building phylogenetic tree ───────────────────────────────────────")

tree_nwk  <- file.path(tmp_dir, "calib_tree.nwk")
ft_bin    <- resolve_bin(FASTTREE_PATH)
cmd_ft    <- sprintf("%s -nt -gtr -quiet %s > %s 2>/dev/null",
                     shQuote(ft_bin), shQuote(aligned_fa), shQuote(tree_nwk))
ft_ret    <- system(cmd_ft)
ft_ok     <- (ft_ret == 0L) && file.exists(tree_nwk) && file.info(tree_nwk)$size > 0L

if (ft_ok) {
  message("  FastTree (GTR) complete")
  tree <- ape::read.tree(tree_nwk)
} else {
  message("  FastTree failed — falling back to NJ (ape)")
  aln_seqs <- read_fasta(aligned_fa)
  mat_list <- lapply(strsplit(toupper(aln_seqs), ""), function(x) x)
  dna_bin  <- ape::as.DNAbin(do.call(rbind, mat_list))
  dm       <- ape::dist.dna(dna_bin, model = "K80", pairwise.deletion = TRUE)
  n_na     <- sum(is.na(as.vector(dm)))
  if (n_na > 0) message("  NOTE: ", n_na, " NA distances — using njs()")
  tree <- ape::njs(dm)
  ape::write.tree(tree, file = tree_nwk)
  message("  NJ tree complete")
}

# Make tip labels safe (ape sometimes mangles them)
tree$tip.label <- gsub("'", "", tree$tip.label)

# Restrict to tips present in taxonomy
common_tips <- intersect(tree$tip.label, rownames(tax_df))
if (length(common_tips) < length(tree$tip.label)) {
  message(sprintf("  Pruning tree: %d → %d tips (keeping only classified)",
                  length(tree$tip.label), length(common_tips)))
  tree <- ape::keep.tip(tree, common_tips)
}

# =============================================================================
# 4. PATRISTIC DISTANCE MATRIX
# =============================================================================
# Distances are sums of branch lengths from the phylogenetic tree.
# When FastTree is used (recommended) they are in units of GTR substitutions
# per site.  When the NJ fallback is used they are K80-corrected distances.
# The two are similar but NOT identical; calibrate and apply using the same
# tree-building method as build_mcrA_db.R (FastTree GTR).
message("\n── 4. Computing patristic distance matrix ──────────────────────────────")

pat_mat <- ape::cophenetic.phylo(tree)
tips    <- rownames(pat_mat)
n_tips  <- length(tips)
message(sprintf("  Matrix: %d × %d (%d unique pairs)",
                n_tips, n_tips, n_tips * (n_tips - 1L) / 2L))

# =============================================================================
# 5. CLASSIFY PAIRS BY DEEPEST DIFFERING RANK
# =============================================================================
message("\n── 5. Classifying sequence pairs by taxonomic rank ────────────────────")

# For each ordered pair (i < j), find the deepest rank where they differ.
# E.g. same genus but different species → boundary_rank = "Species"
#      same family but different genus  → boundary_rank = "Genus"
#      different phylum                 → boundary_rank = "Phylum"

tax_sub <- tax_df[tips, , drop = FALSE]

pair_ranks  <- character(0L)
pair_dists  <- numeric(0L)

# Use vectorised approach over upper triangle
idx_i <- which(upper.tri(matrix(0, n_tips, n_tips)), arr.ind = TRUE)[, 1L]
idx_j <- which(upper.tri(matrix(0, n_tips, n_tips)), arr.ind = TRUE)[, 2L]

message("  Evaluating ", length(idx_i), " pairs …")

# Precompute taxonomy matrix for speed
tax_mat <- as.matrix(tax_sub)   # n_tips × 7

boundary_rank <- vapply(seq_along(idx_i), function(k) {
  ti <- idx_i[k]; tj <- idx_j[k]
  differ <- which(tax_mat[ti, ] != tax_mat[tj, ])
  if (length(differ) == 0L) return("Same")       # identical lineage
  TAX_LEVELS[min(differ)]                         # deepest differing rank
}, character(1L))

pair_dists <- pat_mat[cbind(idx_i, idx_j)]

pairs_df <- data.frame(
  tip_i          = tips[idx_i],
  tip_j          = tips[idx_j],
  dist           = pair_dists,
  boundary_rank  = boundary_rank,
  stringsAsFactors = FALSE
)
# "Same" pairs (identical full lineage) are retained as intra-species negatives
n_same <- sum(pairs_df$boundary_rank == "Same")
message(sprintf("  Informative pairs: %d  (+ %d identical-lineage pairs used for Species)",
                sum(pairs_df$boundary_rank != "Same"), n_same))

# =============================================================================
# 6. YOUDEN-J THRESHOLD PER RANK
# =============================================================================
# For rank R:
#   positive (label = 1) = pairs whose boundary_rank == R  (inter-rank)
#   negative (label = 0) = pairs whose boundary_rank is a rank BELOW R
#                          (within R, i.e. they agree at rank R)
#
# The optimal threshold d* minimises classification error for the question:
#   "does this pair belong to the same R-level taxon?"
# =============================================================================
message("\n── 6. Youden-J threshold estimation ────────────────────────────────────")

rank_order <- c("Species", "Genus", "Family", "Order", "Class", "Phylum")

youden_threshold <- function(dist_vec, label_vec) {
  # label_vec: 1 = inter-rank (positive), 0 = intra-rank (negative)
  # Returns list: threshold, J, sensitivity, specificity
  # NOTE: pROC::coords() "youden" = sensitivity + specificity (raw sum, not -1).
  #       We subtract 1 here to get the standard Youden index J ∈ [0, 1].
  roc_obj <- pROC::roc(response  = label_vec,
                        predictor = dist_vec,
                        direction = "<",   # larger dist → more likely inter-rank
                        quiet     = TRUE)
  coords  <- pROC::coords(roc_obj, x = "best", best.method = "youden",
                           ret = c("threshold", "sensitivity", "specificity",
                                   "youden"))
  # coords may return multiple rows if tied; take first
  row <- as.list(coords[1L, ])
  row$youden <- row$youden - 1   # convert sens+spec → J = sens+spec-1
  row
}

bootstrap_ci <- function(dist_vec, label_vec, n_boot, alpha = 0.05) {
  boot_thresholds <- numeric(n_boot)
  n <- length(dist_vec)
  for (b in seq_len(n_boot)) {
    idx <- sample(n, n, replace = TRUE)
    tryCatch({
      r <- pROC::roc(response  = label_vec[idx],
                     predictor = dist_vec[idx],
                     direction = "<", quiet = TRUE)
      co <- pROC::coords(r, "best", best.method = "youden",
                          ret = "threshold")
      boot_thresholds[b] <- as.numeric(co[1L, "threshold"])
    }, error = function(e) {
      boot_thresholds[b] <<- NA_real_
    })
  }
  boot_thresholds <- boot_thresholds[!is.na(boot_thresholds)]
  quantile(boot_thresholds, probs = c(alpha / 2, 1 - alpha / 2))
}

results <- list()

MAX_IMBALANCE_RATIO <- 10L   # downsample majority class if ratio exceeds this

for (rank in CALIBRATE_RANKS) {
  rank_idx    <- which(rank_order == rank)
  below_ranks <- rank_order[seq_len(rank_idx - 1L)]   # ranks BELOW (finer)

  inter_pairs <- pairs_df$dist[pairs_df$boundary_rank == rank]

  # Species: no finer ranks exist, so use identical-lineage pairs as negatives
  if (rank == "Species") {
    intra_pairs <- pairs_df$dist[pairs_df$boundary_rank == "Same"]
  } else {
    intra_pairs <- pairs_df$dist[pairs_df$boundary_rank %in% below_ranks]
  }

  n_inter <- length(inter_pairs)
  n_intra <- length(intra_pairs)

  if (n_inter < MIN_PAIRS || n_intra < MIN_PAIRS) {
    message(sprintf("  %-8s : SKIPPED (inter=%d, intra=%d; need ≥%d each)",
                    rank, n_inter, n_intra, MIN_PAIRS))
    next
  }

  # Downsample majority class if severely imbalanced (ratio > MAX_IMBALANCE_RATIO)
  # Severe imbalance biases Youden-J toward the majority class.
  ratio <- max(n_inter, n_intra) / min(n_inter, n_intra)
  if (ratio > MAX_IMBALANCE_RATIO) {
    target_n <- min(n_inter, n_intra) * MAX_IMBALANCE_RATIO
    if (n_inter > n_intra) {
      inter_pairs <- sample(inter_pairs, target_n)
      message(sprintf("  %-8s : downsampled inter %d → %d (imbalance ratio %.0f:1)",
                      rank, n_inter, target_n, ratio))
    } else {
      intra_pairs <- sample(intra_pairs, target_n)
      message(sprintf("  %-8s : downsampled intra %d → %d (imbalance ratio %.0f:1)",
                      rank, n_intra, target_n, ratio))
    }
    n_inter <- length(inter_pairs)
    n_intra <- length(intra_pairs)
  }

  dist_vec  <- c(inter_pairs, intra_pairs)
  label_vec <- c(rep(1L, n_inter), rep(0L, n_intra))

  res   <- youden_threshold(dist_vec, label_vec)
  ci    <- bootstrap_ci(dist_vec, label_vec, n_boot = N_BOOT)
  auc   <- as.numeric(pROC::auc(pROC::roc(label_vec, dist_vec,
                                            direction = "<", quiet = TRUE)))

  results[[rank]] <- data.frame(
    rank        = rank,
    n_inter     = n_inter,
    n_intra     = n_intra,
    threshold   = round(res$threshold, 4L),
    ci_lo       = round(ci[1L],        4L),
    ci_hi       = round(ci[2L],        4L),
    sensitivity = round(res$sensitivity,  3L),
    specificity = round(res$specificity,  3L),
    youden_J    = round(res$youden,       3L),
    AUC         = round(auc,              3L),
    stringsAsFactors = FALSE
  )

  message(sprintf("  %-8s : d* = %.4f  [%.4f, %.4f]  J=%.3f  AUC=%.3f  (n_inter=%d, n_intra=%d)",
                  rank,
                  results[[rank]]$threshold,
                  results[[rank]]$ci_lo,
                  results[[rank]]$ci_hi,
                  results[[rank]]$youden_J,
                  results[[rank]]$AUC,
                  n_inter, n_intra))
}

if (length(results) == 0L)
  stop("No ranks had sufficient pairs for calibration. ",
       "Check that tax4mcrA.taxonomy has diverse, fully-classified entries.")

results_df <- do.call(rbind, results)
rownames(results_df) <- NULL

# =============================================================================
# 7. WRITE TSV
# =============================================================================
message("\n── 7. Saving results ────────────────────────────────────────────────────")

write.table(results_df, OUT_TSV, sep = "\t", quote = FALSE, row.names = FALSE)
message("  Thresholds written to: ", OUT_TSV)

# Print ready-to-paste R vector
cat("\n# ── Paste into build_mcrA_db.R ────────────────────────────────────────\n")
cat("# Patristic distances (GTR substitutions/site via FastTree)\n")
cat("PHYLO_RANK_THRESHOLDS <- c(\n")
for (i in seq_len(nrow(results_df))) {
  comma <- if (i < nrow(results_df)) "," else ""
  cat(sprintf("  %-8s = %.4f%s   # J=%.3f  AUC=%.3f  95%%CI [%.4f, %.4f]\n",
              tolower(results_df$rank[i]),
              results_df$threshold[i],
              comma,
              results_df$youden_J[i],
              results_df$AUC[i],
              results_df$ci_lo[i],
              results_df$ci_hi[i]))
}
cat(")\n")

# =============================================================================
# 8. PLOTS
# =============================================================================
message("\n── 8. Generating diagnostic plots ──────────────────────────────────────")

cairo_pdf(OUT_PDF, width = 14, height = 2 * ceiling(length(results) / 2), onefile = TRUE)

for (rank in names(results)) {
  rank_idx    <- which(rank_order == rank)
  below_ranks <- rank_order[seq_len(rank_idx - 1L)]

  inter_d <- pairs_df$dist[pairs_df$boundary_rank == rank]

  # Mirror the same intra-rank logic used in Phase 6:
  # Species has no finer ranks, so "Same" (identical lineage) pairs serve as negatives.
  if (rank == "Species") {
    intra_d <- pairs_df$dist[pairs_df$boundary_rank == "Same"]
  } else {
    intra_d <- pairs_df$dist[pairs_df$boundary_rank %in% below_ranks]
  }

  # Guard: skip plot if either class is empty (shouldn't happen for calibrated ranks,
  # but prevents pROC "response must have two levels" if data are missing).
  if (length(inter_d) == 0L || length(intra_d) == 0L) {
    message(sprintf("  [plot] %-8s : skipped (inter=%d, intra=%d)",
                    rank, length(inter_d), length(intra_d)))
    next
  }

  d_all   <- c(inter_d, intra_d)
  lab_all <- c(rep("inter-rank", length(inter_d)),
                rep("intra-rank", length(intra_d)))
  thr     <- results[[rank]]$threshold

  plot_df <- data.frame(dist = d_all, label = lab_all)

  # ── Density plot ──────────────────────────────────────────────────────────
  p_dens <- ggplot(plot_df, aes(x = dist, fill = label)) +
    geom_density(alpha = 0.5, colour = NA) +
    geom_vline(xintercept = thr, linetype = "dashed", colour = "firebrick",
               linewidth = 0.8) +
    annotate("text", x = thr, y = Inf,
             label = sprintf("d* = %.4f", thr),
             hjust = -0.1, vjust = 1.5, colour = "firebrick", size = 3.5) +
    scale_fill_manual(values = c("inter-rank" = "#E07B54",
                                  "intra-rank" = "#5B8DB8")) +
    labs(title = sprintf("%s threshold (Youden-J)", rank),
         subtitle = sprintf("J = %.3f   AUC = %.3f   n_inter = %d   n_intra = %d",
                             results[[rank]]$youden_J,
                             results[[rank]]$AUC,
                             results[[rank]]$n_inter,
                             results[[rank]]$n_intra),
         x = "Patristic distance (cophenetic)", y = "Density",
         fill = NULL) +
    theme_bw(base_size = 11) +
    theme(legend.position = "top")

  # ── ROC curve ─────────────────────────────────────────────────────────────
  dist_vec  <- c(inter_d, intra_d)
  label_vec <- c(rep(1L, length(inter_d)), rep(0L, length(intra_d)))
  # Guard: pROC requires exactly two levels in response.
  if (length(unique(label_vec)) < 2L) {
    message(sprintf("  [plot] %-8s : ROC skipped (only one label level)", rank))
    next
  }
  roc_obj   <- pROC::roc(label_vec, dist_vec, direction = "<", quiet = TRUE)
  roc_df    <- data.frame(
    fpr = 1 - roc_obj$specificities,
    tpr = roc_obj$sensitivities
  )
  auc_val <- round(as.numeric(pROC::auc(roc_obj)), 3L)

  p_roc <- ggplot(roc_df, aes(x = fpr, y = tpr)) +
    geom_line(colour = "#2C6E49", linewidth = 0.9) +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", colour = "grey50") +
    labs(title = sprintf("%s — ROC curve (AUC = %.3f)", rank, auc_val),
         x = "False positive rate (1 − specificity)",
         y = "True positive rate (sensitivity)") +
    theme_bw(base_size = 11)

  gridExtra::grid.arrange(p_dens, p_roc, ncol = 2L)
}

dev.off()
message("  Plots written to: ", OUT_PDF)

message("\n── Done ─────────────────────────────────────────────────────────────────\n")
