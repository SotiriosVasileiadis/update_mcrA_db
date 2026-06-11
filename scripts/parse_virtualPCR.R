# parse_virtualPCR.R
#
# Parses virtualPCR report.out files produced by run_virtualPCR.sh and
# generates a primer-coverage summary table plus per-taxon coverage tables.
#
# Output format (confirmed from actual report.out):
#   - Blocks separated by exactly 50 underscores  (_{50})
#   - Each block: line 1 = "In silico PCR, primers search for: ..."
#                 line 2 = sequence name (FASTA header without ">")
#                 line 3 = "Target sequence length = N nt"
#   - Amplicon found: line matching "^Amplicon size: <N>bp"
#   - Internal alignment blocks separated by 29 underscores (_{29}) — must NOT
#     be used as block delimiters
#
# FASTA header taxonomy format (semicolon-delimited):
#   Domain;Phylum;Class;Order;Family;Genus;Species;Accession[;]
#
# Run AFTER run_virtualPCR.sh has completed.
# Usage:  Rscript scripts/parse_virtualPCR.R       (from repo root)
#         Rscript parse_virtualPCR.R               (from scripts/)
# ─────────────────────────────────────────────────────────────────────────────

# Resolve the directory that contains this script, so all paths are relative to
# it regardless of where the user calls Rscript from.
#   Rscript scripts/parse_virtualPCR.R  (from repo root)
#   Rscript /abs/path/parse_virtualPCR.R (from anywhere)
#   source("parse_virtualPCR.R")        (from RStudio — falls back to getwd())
BASE_DIR <- tryCatch({
  args <- commandArgs(trailingOnly = FALSE)
  f    <- sub("--file=", "", args[startsWith(args, "--file=")])
  if (length(f) == 1L && nchar(f) > 0L)
    dirname(normalizePath(f, mustWork = FALSE))
  else
    getwd()
}, error = function(e) getwd())

# Results directory: results/primer_coverage/ relative to repo root
RES_DIR <- file.path(BASE_DIR, "..", "results", "primer_coverage")

# Databases: databases/<db>/dada2/<db>.fasta relative to repo root
DATABASES <- list(
  mcrA_ncbi_genome_db = file.path(BASE_DIR, "..", "databases", "mcrA_ncbi_genome_db", "dada2", "mcrA_ncbi_genome_db.fasta"),
  mcrA_ncbi_nt_cur_db = file.path(BASE_DIR, "..", "databases", "mcrA_ncbi_nt_cur_db", "dada2", "mcrA_ncbi_nt_cur_db.fasta"),
  mcrA_ncbi_nt_db     = file.path(BASE_DIR, "..", "databases", "mcrA_ncbi_nt_db",     "dada2", "mcrA_ncbi_nt_db.fasta")
)

TAX_LEVELS <- c("domain", "phylum", "class", "order", "family", "genus")

# ── helpers ───────────────────────────────────────────────────────────────────

count_fasta_seqs <- function(path) {
  lines <- readLines(path, warn = FALSE)
  sum(startsWith(lines, ">"))
}

# Returns all FASTA headers (without leading ">")
get_fasta_headers <- function(fasta_path) {
  lines <- readLines(fasta_path, warn = FALSE)
  sub("^>", "", lines[startsWith(lines, ">")])
}

# Converts a character vector of FASTA header strings into a taxonomy data frame.
# Header format: Domain;Phylum;Class;Order;Family;Genus;Species;Accession[;]
parse_taxonomy_df <- function(seq_names) {
  rows <- lapply(seq_names, function(nm) {
    parts <- trimws(strsplit(nm, ";", fixed = TRUE)[[1]])
    parts <- parts[nchar(parts) > 0]   # drop empty trailing element
    length(parts) <- 6                 # keep Domain–Genus; drop Species & Accession
    setNames(as.list(ifelse(is.na(parts), "Unclassified", parts)), TAX_LEVELS)
  })
  df <- do.call(rbind.data.frame, c(rows, list(stringsAsFactors = FALSE)))
  df$seq_name <- seq_names
  df
}

# For one taxonomy level, compute total / amplified / coverage per unique taxon.
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

# Parses one report.out.  Returns a list with coverage counts AND amplified
# sequence names (needed for taxonomy analysis).
parse_report <- function(report_path, fasta_path) {

  if (!file.exists(report_path)) {
    message("    report.out not found — run run_virtualPCR.sh first.")
    return(NULL)
  }

  raw     <- readLines(report_path, warn = FALSE)
  n_fasta <- count_fasta_seqs(fasta_path)

  # ── Split on 50-underscore block separator only ───────────────────────────
  sep_lines <- which(nchar(raw) == 50L & grepl("^_{50}$", raw))

  if (length(sep_lines) == 0) {
    message("    No block separators found. Check report.out format.")
    return(NULL)
  }

  block_starts <- sep_lines + 1L
  block_ends   <- c(sep_lines[-1] - 1L, length(raw))

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
    n_fasta             = n_fasta,
    n_seqs_with_binding = length(unique(seqs_all)),
    n_amplified         = length(unique(seqs_amp)),
    n_no_binding        = n_fasta - length(unique(seqs_all)),
    amp_sizes           = amp_sizes,
    seqs_amp            = unique(seqs_amp)   # <-- for taxonomy analysis
  )
}

# ── main ──────────────────────────────────────────────────────────────────────

message("=== parse_virtualPCR.R ===")
message("Primers: mlas-mod-F (5'-GGYGGTGTMGGDTTCACMCARTA-3')")
message("         mcrA-rev-R (5'-BGCGTAGTTVGGRTAGT-3', first 7 bp excluded)")
message("Results directory: ", RES_DIR, "\n")

rows        <- list()
all_tax     <- list()    # accumulates taxonomy tables across databases for plotting
tax_dir     <- file.path(RES_DIR, "taxonomy")
dir.create(tax_dir, showWarnings = FALSE, recursive = TRUE)

for (db_name in names(DATABASES)) {
  message("─────────────────────────────────────────────")
  message("Database: ", db_name)

  res <- parse_report(
    file.path(RES_DIR, db_name, "report.out"),
    DATABASES[[db_name]]
  )

  if (is.null(res)) {
    rows[[db_name]] <- data.frame(
      database = db_name, total_seqs = NA, seqs_with_primer_binding = NA,
      amplified = NA, not_amplified = NA, no_binding = NA,
      coverage_pct = NA, amp_min_bp = NA, amp_median_bp = NA, amp_max_bp = NA,
      stringsAsFactors = FALSE)
    next
  }

  bp       <- res$amp_sizes
  coverage <- round(100 * res$n_amplified / res$n_fasta, 2)

  message(sprintf("  Total sequences          : %d", res$n_fasta))
  message(sprintf("  Seqs with primer binding : %d", res$n_seqs_with_binding))
  message(sprintf("  Amplified (both primers) : %d  (%.2f%%)", res$n_amplified, coverage))
  message(sprintf("  Not amplified            : %d", res$n_fasta - res$n_amplified))
  message(sprintf("  No primer binding at all : %d", res$n_no_binding))
  if (length(bp) > 0)
    message(sprintf("  Amplicon size            : min=%d  median=%d  max=%d bp",
                    min(bp), as.integer(median(bp)), max(bp)))

  rows[[db_name]] <- data.frame(
    database                 = db_name,
    total_seqs               = res$n_fasta,
    seqs_with_primer_binding = res$n_seqs_with_binding,
    amplified                = res$n_amplified,
    not_amplified            = res$n_fasta - res$n_amplified,
    no_binding               = res$n_no_binding,
    coverage_pct             = coverage,
    amp_min_bp               = if (length(bp) > 0) min(bp)                  else NA,
    amp_median_bp            = if (length(bp) > 0) as.integer(median(bp))   else NA,
    amp_max_bp               = if (length(bp) > 0) max(bp)                  else NA,
    stringsAsFactors         = FALSE
  )

  # ── per-taxon coverage ───────────────────────────────────────────────────
  message("  Computing per-taxon coverage...")

  all_headers <- get_fasta_headers(DATABASES[[db_name]])
  tax_df      <- parse_taxonomy_df(all_headers)

  tax_tables <- lapply(TAX_LEVELS, function(lvl) {
    tbl <- taxon_coverage(tax_df, res$seqs_amp, lvl)
    tbl$database <- db_name
    tbl[, c("database", "level", "taxon", "total_seqs",
            "amplified", "not_amplified", "coverage_pct")]
  })
  tax_combined <- do.call(rbind, tax_tables)
  rownames(tax_combined) <- NULL

  tax_out <- file.path(tax_dir, paste0("taxonomy_coverage_", db_name, ".tsv"))
  write.table(tax_combined, tax_out, sep = "\t", quote = FALSE, row.names = FALSE)
  message("  Taxonomy coverage saved: ", tax_out)
  all_tax[[db_name]] <- tax_combined

  # Print phylum-level summary to console
  message("\n  ── Phylum-level coverage ──────────────────────")
  phylum_tbl <- tax_combined[tax_combined$level == "phylum", ]
  phylum_tbl <- phylum_tbl[order(-phylum_tbl$coverage_pct, -phylum_tbl$total_seqs), ]
  print(phylum_tbl[, c("taxon","total_seqs","amplified","coverage_pct")],
        row.names = FALSE)

  # Print genus-level summary (top 20 by total_seqs)
  genus_tbl <- tax_combined[tax_combined$level == "genus", ]
  genus_tbl <- genus_tbl[order(-genus_tbl$total_seqs), ]
  message("\n  ── Genus-level coverage (top 20 by total_seqs) ──")
  print(head(genus_tbl[, c("taxon","total_seqs","amplified","coverage_pct")], 20),
        row.names = FALSE)
  message("")
}

message("\n═══════════════════════════════════════════════════════════")
message("COVERAGE SUMMARY  —  mlas-mod-F / mcrA-rev-R")
message("═══════════════════════════════════════════════════════════")
summary_df <- do.call(rbind, rows)
rownames(summary_df) <- NULL
print(summary_df, row.names = FALSE)

out_tsv <- file.path(RES_DIR, "coverage_summary.tsv")
write.table(summary_df, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
message("\nSummary saved      : ", out_tsv)
message("Taxonomy tables in : ", tax_dir)

# ── line plots ────────────────────────────────────────────────────────────────
# y-axis: % of sequences per taxon that were amplified (coverage_pct)
# x-axis: database (categorical, ordered as in DATABASES list)
# one line per taxon; domain excluded; genus/family capped at top 20 by total_seqs

make_single_run_plots <- function(all_tax, db_names, plot_dir,
                                  max_taxa = c(phylum=Inf, class=Inf, order=Inf,
                                               family=20, genus=20)) {
  dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
  plot_levels <- setdiff(TAX_LEVELS, "domain")

  # Combine taxonomy tables from all databases
  # (the 'database' column is already set in each entry by the main loop)
  if (length(all_tax) == 0L) {
    message("  No taxonomy data available — skipping plots.")
    return(invisible(NULL))
  }
  tax_df <- do.call(rbind, all_tax)
  rownames(tax_df) <- NULL
  if (!is.data.frame(tax_df) || nrow(tax_df) == 0L) {
    message("  Empty taxonomy data — skipping plots.")
    return(invisible(NULL))
  }

  db_short  <- sub("mcrA_ncbi_", "", db_names)
  x_vals    <- seq_along(db_names)

  for (lvl in plot_levels) {

    sub <- tax_df[!is.na(tax_df$level) & tax_df$level == lvl, , drop = FALSE]
    if (!is.data.frame(sub) || nrow(sub) == 0L) next

    # Keep top taxa by mean total_seqs across databases
    mt   <- max_taxa[lvl]
    avg  <- tapply(sub$total_seqs, sub$taxon, mean)
    keep <- names(sort(avg, decreasing = TRUE))[seq_len(min(mt, length(avg)))]
    sub  <- sub[sub$taxon %in% keep, ]

    taxa <- unique(sub$taxon)
    n_t  <- length(taxa)

    # Coverage matrix: rows = taxa, cols = databases
    mat <- matrix(NA_real_, n_t, length(db_names),
                  dimnames = list(taxa, db_names))
    for (tx in taxa) {
      r <- sub[sub$taxon == tx, ]
      for (j in seq_along(db_names)) {
        v <- r$coverage_pct[r$database == db_names[j]]
        if (length(v) == 1L) mat[tx, j] <- v
      }
    }

    pal  <- tryCatch(hcl.colors(max(n_t, 2L), "Dark 3"),
                     error = function(e) rainbow(max(n_t, 2L)))
    cols <- pal[seq_len(n_t)]

    pdf_path <- file.path(plot_dir, paste0("coverage_", lvl, ".pdf"))
    pdf(pdf_path, width = 6.5, height = 5.5)
    par(mar = c(6.5, 4.5, 3, 0.5))

    matplot(x_vals, t(mat),
            type = "b", lty = 1, lwd = 1.8, pch = 19, cex = 0.7,
            col  = cols,
            xlab = "", ylab = "Amplified sequences (%)",
            main = paste0("Primer coverage — ", lvl, " level"),
            ylim = c(0, 100), xaxt = "n")
    axis(1, at = x_vals, labels = db_short, las = 2, cex.axis = 0.85)
    grid(nx = NA, ny = NULL, col = "grey88", lty = 1)

    cex_leg  <- max(0.45, min(0.75, 9 / n_t))
    ncol_leg <- if (n_t > 15L) 2L else 1L
    legend("bottomright", legend = taxa,
           col = cols, lty = 1, lwd = 1.5, pch = 19, pt.cex = 0.6,
           cex = cex_leg, bty = "n", ncol = ncol_leg)

    dev.off()
    message("  Plot: ", pdf_path)
  }
}

plot_dir <- file.path(tax_dir, "plots")
message("\nGenerating line plots...")
make_single_run_plots(all_tax, names(DATABASES), plot_dir)
message("Plots in: ", plot_dir)
