# parse_mismatch_sweep.R
#
# Parses all results/primer_coverage_sweep_errorsN/ folders produced by
# run_mismatch_sweep.sh and produces:
#   1. Cross-mismatch coverage comparison (long + wide TSV)
#   2. Per-taxon coverage tables at each mismatch level (one TSV per tax level)
#   3. Multi-panel line plots (single PDF)
#
# FASTA header taxonomy format (semicolon-delimited):
#   Domain;Phylum;Class;Order;Family;Genus;Species;Accession[;]
#
# Usage (from any directory):
#   Rscript scripts/parse_mismatch_sweep.R   (from repo root)
#   Rscript parse_mismatch_sweep.R           (from scripts/)
# ─────────────────────────────────────────────────────────────────────────────

# Resolve the directory that contains this script, so all paths are relative to
# it regardless of where the user calls Rscript from.
#   Rscript scripts/parse_mismatch_sweep.R          (from repo root)
#   Rscript /abs/path/parse_mismatch_sweep.R (from anywhere)
#   source("parse_mismatch_sweep.R")        (from RStudio — falls back to getwd())
BASE_DIR <- tryCatch({
  args <- commandArgs(trailingOnly = FALSE)
  f    <- sub("--file=", "", args[startsWith(args, "--file=")])
  if (length(f) == 1L && nchar(f) > 0L)
    dirname(normalizePath(f, mustWork = FALSE))
  else
    getwd()
}, error = function(e) getwd())

# Repo root (one level up from scripts/)
REPO_DIR <- file.path(BASE_DIR, "..")

# Databases: databases/<db>/dada2/<db>.fasta relative to repo root
DATABASES <- list(
  mcrA_ncbi_genome_db = file.path(REPO_DIR, "databases", "mcrA_ncbi_genome_db", "dada2", "mcrA_ncbi_genome_db.fasta"),
  mcrA_ncbi_nt_cur_db = file.path(REPO_DIR, "databases", "mcrA_ncbi_nt_cur_db", "dada2", "mcrA_ncbi_nt_cur_db.fasta"),
  mcrA_ncbi_nt_db     = file.path(REPO_DIR, "databases", "mcrA_ncbi_nt_db",     "dada2", "mcrA_ncbi_nt_db.fasta")
)

TAX_LEVELS <- c("domain", "phylum", "class", "order", "family", "genus")

# ── helpers (shared with parse_virtualPCR.R) ──────────────────────────────────

count_fasta_seqs <- function(path) {
  sum(startsWith(readLines(path, warn = FALSE), ">"))
}

get_fasta_headers <- function(fasta_path) {
  lines <- readLines(fasta_path, warn = FALSE)
  sub("^>", "", lines[startsWith(lines, ">")])
}

# Converts FASTA header strings to a taxonomy data frame (Domain through Genus).
parse_taxonomy_df <- function(seq_names) {
  rows <- lapply(seq_names, function(nm) {
    parts <- trimws(strsplit(nm, ";", fixed = TRUE)[[1]])
    parts <- parts[nchar(parts) > 0]
    length(parts) <- 6
    setNames(as.list(ifelse(is.na(parts), "Unclassified", parts)), TAX_LEVELS)
  })
  df <- do.call(rbind.data.frame, c(rows, list(stringsAsFactors = FALSE)))
  df$seq_name <- seq_names
  df
}

# Coverage per unique taxon at one taxonomy level.
taxon_coverage <- function(tax_df, amplified_seqs, level) {
  tax_df$is_amp <- tax_df$seq_name %in% amplified_seqs
  taxa  <- tax_df[[level]]
  total <- tapply(rep(1L, nrow(tax_df)), taxa, sum)
  amp   <- tapply(tax_df$is_amp,          taxa, sum)
  data.frame(
    level        = level,
    taxon        = names(total),
    total_seqs   = as.integer(total),
    amplified    = as.integer(amp),
    not_amplified= as.integer(total) - as.integer(amp),
    coverage_pct = round(100 * as.integer(amp) / as.integer(total), 2),
    stringsAsFactors = FALSE,
    row.names    = NULL
  )
}

# Parses one report.out. Returns coverage counts AND amplified sequence names.
parse_report <- function(report_path, fasta_path) {
  if (!file.exists(report_path)) return(NULL)

  raw     <- readLines(report_path, warn = FALSE)
  n_fasta <- count_fasta_seqs(fasta_path)
  sep_idx <- which(nchar(raw) == 50L & grepl("^_{50}$", raw))
  if (length(sep_idx) == 0) return(NULL)

  block_starts <- sep_idx + 1L
  block_ends   <- c(sep_idx[-1] - 1L, length(raw))

  seqs_all  <- character(0)
  seqs_amp  <- character(0)
  amp_sizes <- integer(0)

  for (i in seq_along(block_starts)) {
    block <- raw[block_starts[i]:block_ends[i]]
    block <- block[nchar(trimws(block)) > 0]
    if (length(block) < 2) next
    seq_name <- trimws(block[2])
    seqs_all <- c(seqs_all, seq_name)
    amp_lines <- grep("^Amplicon size:", block, value = TRUE)
    if (length(amp_lines) > 0) {
      seqs_amp  <- c(seqs_amp, seq_name)
      sizes     <- as.integer(regmatches(amp_lines, regexpr("\\d+", amp_lines)))
      amp_sizes <- c(amp_sizes, sizes[!is.na(sizes)])
    }
  }

  list(
    n_fasta     = n_fasta,
    n_binding   = length(unique(seqs_all)),
    n_amplified = length(unique(seqs_amp)),
    amp_sizes   = amp_sizes,
    seqs_amp    = unique(seqs_amp)   # <-- for taxonomy analysis
  )
}

# ── discover results folders ──────────────────────────────────────────────────
# Sweep results live in results/primer_coverage_sweep_errorsN/ under repo root

results_root <- file.path(REPO_DIR, "results")
res_dirs_full <- sort(list.dirs(results_root, recursive = FALSE, full.names = FALSE))
res_dirs <- res_dirs_full[grepl("^primer_coverage_sweep_errors\\d+$", res_dirs_full)]

if (length(res_dirs) == 0)
  stop("No primer_coverage_sweep_errorsN/ folders found under results/. ",
       "Run run_mismatch_sweep.sh first.")

error_levels <- as.integer(sub("primer_coverage_sweep_errors", "", res_dirs))
message("Found mismatch levels: ", paste(error_levels, collapse = ", "))

# Pre-load FASTA taxonomy tables (one per database; same across all mismatch levels)
message("Loading FASTA taxonomy maps...")
fasta_tax <- lapply(names(DATABASES), function(db) {
  parse_taxonomy_df(get_fasta_headers(DATABASES[[db]]))
})
names(fasta_tax) <- names(DATABASES)

# ── parse all results ─────────────────────────────────────────────────────────

all_rows    <- list()   # for overall coverage summary
tax_records <- list()   # for per-taxon analysis: list of data frames

for (k in seq_along(res_dirs)) {
  nerr    <- error_levels[k]
  res_dir <- file.path(results_root, res_dirs[k])
  message("\n── number3errors = ", nerr, "  (", res_dirs[k], ") ──")

  for (db_name in names(DATABASES)) {
    report <- file.path(res_dir, db_name, "report.out")
    res    <- parse_report(report, DATABASES[[db_name]])

    if (is.null(res)) {
      message("  ", db_name, ": report.out not found — skipping")
      next
    }

    bp       <- res$amp_sizes
    coverage <- round(100 * res$n_amplified / res$n_fasta, 2)
    message(sprintf("  %-25s  amplified=%d / %d  (%.2f%%)",
                    db_name, res$n_amplified, res$n_fasta, coverage))

    all_rows[[paste(nerr, db_name, sep = "_")]] <- data.frame(
      number3errors      = nerr,
      database           = db_name,
      total_seqs         = res$n_fasta,
      seqs_with_binding  = res$n_binding,
      amplified          = res$n_amplified,
      not_amplified      = res$n_fasta - res$n_amplified,
      coverage_pct       = coverage,
      amp_min_bp         = if (length(bp) > 0) min(bp)                 else NA_integer_,
      amp_median_bp      = if (length(bp) > 0) as.integer(median(bp))  else NA_integer_,
      amp_max_bp         = if (length(bp) > 0) max(bp)                 else NA_integer_,
      stringsAsFactors   = FALSE
    )

    # ── per-taxon rows for this (nerr, db_name) ─────────────────────────────
    tax_df <- fasta_tax[[db_name]]
    for (lvl in TAX_LEVELS) {
      tbl <- taxon_coverage(tax_df, res$seqs_amp, lvl)
      tbl$database      <- db_name
      tbl$number3errors <- nerr
      tax_records[[paste(nerr, db_name, lvl, sep = "_")]] <- tbl
    }
  }
}

# ── overall coverage summary ──────────────────────────────────────────────────

summary_df <- do.call(rbind, all_rows)
rownames(summary_df) <- NULL
summary_df <- summary_df[order(summary_df$database, summary_df$number3errors), ]

message("\n═══════════════════════════════════════════════════════════")
message("MISMATCH SWEEP — mlas-mod-F / mcrA-rev-R")
message("  rev trimmed: 5'-BGCGTAGTTVGGRTAGT-3' (7 bp excluded at 5' end)")
message("═══════════════════════════════════════════════════════════")
print(summary_df[, c("number3errors","database","total_seqs",
                     "amplified","coverage_pct",
                     "amp_min_bp","amp_median_bp","amp_max_bp")],
      row.names = FALSE)

wide <- reshape(
  summary_df[, c("database","number3errors","amplified","coverage_pct")],
  idvar     = "database",
  timevar   = "number3errors",
  v.names   = c("amplified","coverage_pct"),
  direction = "wide"
)
wide <- wide[order(wide$database), ]

message("\n── Coverage (%) by mismatch level ──────────────────────────")
print(wide, row.names = FALSE)

# ── per-taxon coverage tables ─────────────────────────────────────────────────

tax_all <- do.call(rbind, tax_records)
rownames(tax_all) <- NULL

# Reorder columns for readability
tax_all <- tax_all[, c("number3errors", "database", "level", "taxon",
                        "total_seqs", "amplified", "not_amplified", "coverage_pct")]
tax_all <- tax_all[order(tax_all$level, tax_all$database,
                          tax_all$number3errors, -tax_all$total_seqs), ]

message("\n═══════════════════════════════════════════════════════════")
message("PER-TAXON COVERAGE")
message("═══════════════════════════════════════════════════════════")

# Print phylum-level table to console (most informative at a glance)
message("\n── Phylum-level coverage across mismatch levels ───────────")
phylum_wide <- reshape(
  tax_all[tax_all$level == "phylum",
          c("database","number3errors","taxon","coverage_pct")],
  idvar     = c("database","taxon"),
  timevar   = "number3errors",
  v.names   = "coverage_pct",
  direction = "wide"
)
phylum_wide <- phylum_wide[order(phylum_wide$database, -rowSums(phylum_wide[, -(1:2), drop=FALSE], na.rm=TRUE)), ]
print(phylum_wide, row.names = FALSE)

# Print genus-level for each db (top 15 by total_seqs at max mismatch level)
for (db in names(DATABASES)) {
  max_nerr   <- max(error_levels)
  genus_tbl  <- tax_all[tax_all$level == "genus" &
                         tax_all$database == db &
                         tax_all$number3errors == max_nerr, ]
  genus_tbl  <- genus_tbl[order(-genus_tbl$total_seqs), ]
  message(sprintf("\n── Genus-level, %s, errors=%d (top 15 by total_seqs) ──", db, max_nerr))
  print(head(genus_tbl[, c("taxon","total_seqs","amplified","coverage_pct")], 15),
        row.names = FALSE)
}

# ── save ──────────────────────────────────────────────────────────────────────

out_long <- file.path(results_root, "mismatch_sweep_summary.tsv")
out_wide <- file.path(results_root, "mismatch_sweep_wide.tsv")
write.table(summary_df, out_long, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(wide,       out_wide, sep = "\t", quote = FALSE, row.names = FALSE)

# Save one TSV per taxonomy level (all databases × all mismatch levels combined)
tax_dir <- file.path(results_root, "taxonomy")
dir.create(tax_dir, showWarnings = FALSE, recursive = TRUE)

for (lvl in TAX_LEVELS) {
  out_tax <- file.path(tax_dir, paste0("sweep_taxonomy_", lvl, ".tsv"))
  write.table(
    tax_all[tax_all$level == lvl, ],
    out_tax, sep = "\t", quote = FALSE, row.names = FALSE
  )
}

out_tax_all <- file.path(tax_dir, "sweep_taxonomy_all_levels.tsv")
tax_out <- tax_all[tax_all$level != "domain", ]
# rename columns
names(tax_out)[names(tax_out) == "number3errors"] <- "3-end mismatches"
names(tax_out)[names(tax_out) == "amplified"]     <- "hits"
names(tax_out)[names(tax_out) == "not_amplified"] <- "missed"
# rename database values
tax_out$database <- c(
  mcrA_ncbi_genome_db = "MCRA_FULL",
  mcrA_ncbi_nt_cur_db = "mlas-MCRA_CUR",
  mcrA_ncbi_nt_db     = "mlas-MCRA"
)[tax_out$database]
write.table(tax_out, out_tax_all, sep = "\t", quote = FALSE, row.names = FALSE)

message("\nLong table         : ", out_long)
message("Wide table         : ", out_wide)
message("Taxonomy tables in : ", tax_dir)
message("  Files: ", paste0("sweep_taxonomy_", TAX_LEVELS, ".tsv", collapse = ", "))
message("  Combined: ", out_tax_all)

# ── line plots ────────────────────────────────────────────────────────────────
# y-axis: % of sequences per taxon that were amplified (coverage_pct)
# x-axis: number3errors; one line per taxon; one panel per database
# Domain is excluded. Genus/family capped at top 20 taxa by mean total_seqs.

make_sweep_plots <- function(tax_all, db_names, plot_dir,
                             max_taxa = c(phylum=Inf, class=Inf, order=Inf,
                                          family=20, genus=20)) {
  dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
  plot_levels <- setdiff(TAX_LEVELS, "domain")
  n_db  <- length(db_names)
  n_lvl <- length(plot_levels)

  db_labels <- c(mcrA_ncbi_genome_db = "MCRA_FULL",
                 mcrA_ncbi_nt_cur_db = "mlas-MCRA_CUR",
                 mcrA_ncbi_nt_db     = "mlas-MCRA")

  # All levels x databases on a single page, plus a shared legend column per row.
  # Layout: n_lvl rows x (n_db + 1) cols; last column is narrower (legend only).
  pdf_path <- file.path(plot_dir, "sweep_taxonomy_all_levels.pdf")
  pdf(pdf_path,
      width  = 4.5 * n_db + 3.5,
      height = 4.2 * n_lvl + 0.6)

  par(oma = c(1, 0, 2, 0.5))
  lay_mat <- matrix(seq_len(n_lvl * (n_db + 1L)), nrow = n_lvl, byrow = TRUE)
  layout(lay_mat, widths = c(rep(4, n_db), 2.5))

  data_mar   <- c(5.5, 8, 2, 0.5)
  legend_mar <- c(0.5, 0.5, 0.5, 0.5)

  for (lvl in plot_levels) {

    # ── global colour map for this level (consistent across all databases) ────
    all_sub <- tax_all[tax_all$level == lvl, ]
    mt      <- max_taxa[lvl]
    avg     <- tapply(all_sub$total_seqs, all_sub$taxon, mean)
    keep    <- names(sort(avg, decreasing = TRUE))[seq_len(min(mt, length(avg)))]
    n_t_all <- length(keep)
    pal      <- tryCatch(hcl.colors(max(n_t_all, 2L), "Dark 3"),
                         error = function(e) rainbow(max(n_t_all, 2L)))
    cols_all <- setNames(pal[seq_len(n_t_all)], keep)

    # ── database panels ───────────────────────────────────────────────────────
    for (db in db_names) {
      par(mar = data_mar, cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5)

      sub <- tax_all[tax_all$level == lvl & tax_all$database == db, ]
      sub <- sub[sub$taxon %in% keep, ]

      if (nrow(sub) == 0L) { plot.new(); next }

      # Taxa present in this db, preserved in global order for colour consistency
      taxa   <- keep[keep %in% unique(sub$taxon)]
      n_t    <- length(taxa)
      x_vals <- sort(unique(sub$number3errors))

      mat <- matrix(NA_real_, n_t, length(x_vals),
                    dimnames = list(taxa, as.character(x_vals)))
      for (tx in taxa) {
        r <- sub[sub$taxon == tx, ]
        for (j in seq_along(x_vals)) {
          v <- r$coverage_pct[r$number3errors == x_vals[j]]
          if (length(v) == 1L) mat[tx, j] <- v
        }
      }

      cols        <- cols_all[taxa]
      panel_title <- if (lvl == plot_levels[1]) db_labels[db] else ""
      y_label     <- if (db == db_names[1]) paste0(lvl, "\nhits (%)") else ""

      matplot(x_vals, t(mat),
              type = "b", lty = 1, lwd = 1.8, pch = 19, cex = 0.98,
              col  = cols,
              xlab = if (lvl == plot_levels[n_lvl]) "3'-end allowed mismatches" else "",
              ylab = y_label,
              main = panel_title,
              ylim = c(0, 100), xaxt = "n")
      axis(1, at = x_vals)
      grid(nx = NA, ny = NULL, col = "grey88", lty = 1)
    }

    # ── shared legend panel for this level ────────────────────────────────────
    par(mar = legend_mar)
    plot.new()
    cex_leg  <- 1.08
    ncol_leg <- 1L
    legend("center", legend = keep,
           col = cols_all[keep], lty = 1, lwd = 1.5, pch = 19, pt.cex = 0.83,
           cex = cex_leg, bty = "n", ncol = ncol_leg)
  }

  mtext("Primer coverage by mismatch tolerance  —  mlas-mod-F / mcrA-rev-R",
        outer = TRUE, side = 3, cex = 1.5, font = 2)
  dev.off()
  message("  Plot: ", pdf_path)
}

plot_dir <- file.path(tax_dir, "plots")
message("\nGenerating line plots...")
db_plot_order <- c("mcrA_ncbi_nt_db", "mcrA_ncbi_nt_cur_db", "mcrA_ncbi_genome_db")
make_sweep_plots(tax_all, db_plot_order, plot_dir)
message("Plots in: ", plot_dir)
