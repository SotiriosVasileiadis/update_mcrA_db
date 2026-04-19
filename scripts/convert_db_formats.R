#!/usr/bin/env Rscript
# =============================================================================
# convert_db_formats.R
# =============================================================================
# Converts mcrA DADA2-format reference databases into QIIME2 and Mothur formats.
#
# DADA2 header format (input):
#   >Domain;Phylum;Class;Order;Family;Genus;Species;Accession[;]
#
# Output folder structure (created inside the databases directory):
#   <db_basename>/
#     dada2/   <db_basename>.fasta           (copy of original)
#     qiime2/  <db_basename>_seqs.fasta      (accession-only headers)
#              <db_basename>_taxonomy.tsv    (Feature ID <tab> k__...;p__...;...)
#     mothur/  <db_basename>.fasta           (accession-only headers)
#              <db_basename>.taxonomy        (SeqID <tab> Domain;...;Species;)
#
# Usage (from terminal):
#   Rscript convert_db_formats.R /path/to/databases
#
# If no argument is given, the script uses its own directory.
# =============================================================================

QIIME2_PREFIXES <- c("k__", "p__", "c__", "o__", "f__", "g__", "s__")

# -----------------------------------------------------------------------------
# Read a FASTA file into a data.frame with columns: header, sequence
# -----------------------------------------------------------------------------
read_fasta <- function(filepath) {
  lines <- readLines(filepath, warn = FALSE)
  is_header <- startsWith(lines, ">")
  header_idx <- which(is_header)

  n <- length(header_idx)
  headers   <- character(n)
  sequences <- character(n)

  for (i in seq_len(n)) {
    start <- header_idx[i] + 1L
    end   <- if (i < n) header_idx[i + 1L] - 1L else length(lines)
    headers[i]   <- lines[header_idx[i]]
    # Guard against start > end (empty sequence block) – R's 1:0 is NOT empty
    seq_block    <- if (start <= end) lines[start:end] else character(0L)
    sequences[i] <- paste(seq_block, collapse = "")
  }

  data.frame(header = headers, sequence = sequences, stringsAsFactors = FALSE)
}

# -----------------------------------------------------------------------------
# Parse a DADA2 header into accession + 7-element taxonomy vector
# Returns a list(accession, taxonomy) or stops with an informative message
# -----------------------------------------------------------------------------
parse_dada2_header <- function(header) {
  h     <- sub("^>", "", header)      # strip leading ">"
  h     <- sub(";$", "", h)           # strip optional trailing ";"
  parts <- strsplit(h, ";")[[1]]

  if (length(parts) < 8L) {
    stop(sprintf(
      "Expected >=8 semicolon-separated fields, got %d: %s",
      length(parts), header
    ))
  }

  list(
    accession = trimws(parts[8L]),
    taxonomy  = parts[1:7]            # Domain … Species
  )
}

# -----------------------------------------------------------------------------
# Make accession IDs unique by appending _1, _2, ... to duplicates
# Singletons keep their original ID unchanged
# -----------------------------------------------------------------------------
make_unique_ids <- function(accessions) {
  counts <- table(accessions)
  seen   <- integer(length(counts))
  names(seen) <- names(counts)

  ids <- character(length(accessions))
  for (i in seq_along(accessions)) {
    acc <- accessions[i]
    if (counts[acc] == 1L) {
      ids[i] <- acc
    } else {
      seen[acc] <- seen[acc] + 1L
      ids[i] <- paste0(acc, "_", seen[acc])
    }
  }
  ids
}

# -----------------------------------------------------------------------------
# Convert a single DADA2 FASTA → DADA2 / QIIME2 / Mothur under out_dir
# -----------------------------------------------------------------------------
convert_database <- function(fasta_path, out_dir) {
  basename_db <- tools::file_path_sans_ext(basename(fasta_path))

  # ── read ──────────────────────────────────────────────────────────────────
  fa <- read_fasta(fasta_path)

  # ── parse headers ─────────────────────────────────────────────────────────
  parsed   <- vector("list", nrow(fa))
  bad_rows <- integer(0)

  for (i in seq_len(nrow(fa))) {
    result <- tryCatch(
      parse_dada2_header(fa$header[i]),
      error = function(e) {
        message("  WARNING – row ", i, ": ", conditionMessage(e))
        NULL
      }
    )
    if (is.null(result)) {
      bad_rows <- c(bad_rows, i)
    } else {
      parsed[[i]] <- result
    }
  }

  if (length(bad_rows) > 0L) {
    fa      <- fa[-bad_rows, ]
    parsed  <- parsed[-bad_rows]
    message("  Skipped ", length(bad_rows), " unparseable header(s).")
  }

  accessions <- vapply(parsed, `[[`, character(1L), "accession")
  taxonomies <- do.call(rbind, lapply(parsed, `[[`, "taxonomy"))  # n × 7 matrix

  # ── deduplicate IDs ───────────────────────────────────────────────────────
  unique_ids <- make_unique_ids(accessions)

  # ── create subdirectories ─────────────────────────────────────────────────
  dada2_dir  <- file.path(out_dir, "dada2")
  qiime2_dir <- file.path(out_dir, "qiime2")
  mothur_dir <- file.path(out_dir, "mothur")
  for (d in c(dada2_dir, qiime2_dir, mothur_dir)) {
    dir.create(d, showWarnings = FALSE, recursive = TRUE)
  }

  # ── DADA2 – copy original ─────────────────────────────────────────────────
  file.copy(fasta_path, file.path(dada2_dir, basename(fasta_path)),
            overwrite = TRUE)

  # ── QIIME2 ────────────────────────────────────────────────────────────────
  #   seqs.fasta  : >UniqueID
  #   taxonomy.tsv: Feature ID \t k__D;p__P;c__C;o__O;f__F;g__G;s__S
  qiime2_seqs <- file.path(qiime2_dir, paste0(basename_db, "_seqs.fasta"))
  qiime2_tax  <- file.path(qiime2_dir, paste0(basename_db, "_taxonomy.tsv"))

  fasta_lines <- character(nrow(fa) * 2L)
  fasta_lines[seq(1, length(fasta_lines), 2)] <- paste0(">", unique_ids)
  fasta_lines[seq(2, length(fasta_lines), 2)] <- fa$sequence
  writeLines(fasta_lines, qiime2_seqs)

  taxon_strings <- apply(taxonomies, 1, function(tax) {
    paste(paste0(QIIME2_PREFIXES, tax), collapse = ";")
  })
  tax_df <- data.frame(
    `Feature ID` = unique_ids,
    Taxon        = taxon_strings,
    check.names  = FALSE,
    stringsAsFactors = FALSE
  )
  write.table(tax_df, qiime2_tax,
              sep = "\t", quote = FALSE, row.names = FALSE)

  # ── Mothur ────────────────────────────────────────────────────────────────
  #   seqs.fasta  : >UniqueID
  #   taxonomy    : UniqueID \t Domain;Phylum;Class;Order;Family;Genus;Species;
  mothur_seqs <- file.path(mothur_dir, paste0(basename_db, ".fasta"))
  mothur_tax  <- file.path(mothur_dir, paste0(basename_db, ".taxonomy"))

  writeLines(fasta_lines, mothur_seqs)   # same FASTA as QIIME2

  # Spaces → underscores in taxon names (required by Mothur)
  tax_clean <- apply(taxonomies, 2, function(col) gsub(" ", "_", col))
  mothur_taxon_strings <- apply(tax_clean, 1, function(tax) {
    paste0(paste(tax, collapse = ";"), ";")   # trailing ";" required
  })
  mothur_lines <- paste(unique_ids, mothur_taxon_strings, sep = "\t")
  writeLines(mothur_lines, mothur_tax)

  nrow(fa)
}

# =============================================================================
# MAIN
# =============================================================================
args   <- commandArgs(trailingOnly = TRUE)
db_dir <- if (length(args) >= 1L) args[1L] else dirname(normalizePath(sys.frame(0)$ofile, mustWork = FALSE))

# Fallback when called interactively or path resolution fails
if (!nzchar(db_dir) || !dir.exists(db_dir)) {
  db_dir <- getwd()
}

fasta_files <- sort(list.files(db_dir, pattern = "\\.fasta$",
                               full.names = FALSE, recursive = FALSE))

if (length(fasta_files) == 0L) {
  stop("No .fasta files found in: ", db_dir)
}

cat(sprintf("Found %d database(s) in: %s\n\n", length(fasta_files), db_dir))

for (fname in fasta_files) {
  fasta_path  <- file.path(db_dir, fname)
  basename_db <- tools::file_path_sans_ext(fname)
  out_dir     <- file.path(db_dir, basename_db)

  cat(sprintf("Processing: %s\n", fname))
  n <- convert_database(fasta_path, out_dir)
  cat(sprintf("  Sequences converted : %s\n", format(n, big.mark = ",")))
  cat(sprintf("  Output folder       : %s/\n", out_dir))
  cat(sprintf("    dada2/  %s.fasta\n",              basename_db))
  cat(sprintf("    qiime2/ %s_seqs.fasta\n",         basename_db))
  cat(sprintf("    qiime2/ %s_taxonomy.tsv\n",       basename_db))
  cat(sprintf("    mothur/ %s.fasta\n",              basename_db))
  cat(sprintf("    mothur/ %s.taxonomy\n\n",         basename_db))
}

cat("Done.\n")
