#!/usr/bin/env Rscript
# =============================================================================
# build_mcrA_db_ncbi_annot.R
#
# Builds a DADA2-formatted mcrA reference database directly from NCBI genome
# annotations, without relying on BLAST against a seed template.
#
# Workflow:
#   1. Query NCBI Datasets v2 API for all methanogenic archaeal genome
#      accessions across the major methanogenic lineages
#   2. Fetch complete taxonomic lineages per taxon via NCBI Taxonomy (rentrez)
#   3. Download annotated CDS nucleotide FASTA per genome (datasets CLI / FTP)
#   4. Extract mcrA-annotated CDS by matching FASTA header patterns:
#        [gene=mcrA...]  OR  methyl-coenzyme M reductase * alpha
#   5. Translate each CDS and run hmmsearch against TIGR03256.hmm (protein HMM)
#   6. Retain sequences covering >= HMM_MIN_COVERAGE (default 90 %) of profile
#   7. Trim nucleotide sequences to the HMM envelope coordinates (nt)
#   8. Deduplicate by exact trimmed nucleotide sequence
#   9. Format DADA2 taxonomy header and write final database FASTA
#
# Output:
#   mcrA_ncbi_genome_db.fasta  — DADA2-format database (header = taxonomy path)
#   mcrA_ncbi_genome_stats.tsv — per-rank sequence counts
#
# Dependencies (R):   httr, jsonlite, rentrez, Biostrings, dplyr, stringr
# Dependencies (CLI): hmmer (hmmsearch + hmmpress), ncbi-datasets-cli (optional)
#
# Usage:
#   Rscript build_mcrA_db_ncbi_annot.R
#   (run from the folder containing TIGR03256.hmm)
# =============================================================================

# ── 0.  Dependencies ──────────────────────────────────────────────────────────
required_pkgs <- c("httr", "jsonlite", "rentrez", "Biostrings", "dplyr", "stringr")
missing_pkgs  <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0L) {
  stop(
    "Install missing packages first:\n",
    "  BiocManager::install('Biostrings')\n",
    "  install.packages(c(",
    paste(sprintf('"%s"', setdiff(missing_pkgs, "Biostrings")), collapse = ", "), "))"
  )
}

suppressPackageStartupMessages({
  library(httr);      library(jsonlite);  library(rentrez)
  library(Biostrings);library(dplyr);     library(stringr)
})

# ── 1.  Configuration ─────────────────────────────────────────────────────────
HMM_PROFILE      <- "TIGR03256.hmm"          # must exist in working directory

OUT_FASTA        <- "mcrA_ncbi_genome_db.fasta"
OUT_STATS        <- "mcrA_ncbi_genome_stats.tsv"
CACHE_DIR        <- "ncbi_annot_cache"
GENOME_CDS_DIR   <- file.path(CACHE_DIR, "cds")
MCRA_RAW_DIR     <- file.path(CACHE_DIR, "mcra_raw")   # per-genome extracted mcrA nt
ACCESSION_CACHE  <- file.path(CACHE_DIR, "genome_accessions.rds")
ASM_FTP_CACHE    <- file.path(CACHE_DIR, "assembly_ftp_map.rds")

# Completeness threshold: fraction of TIGR03256 profile that must be covered
HMM_MIN_COVERAGE <- 0.90

# Minimum nucleotide length of a trimmed mcrA to include in the database
MIN_NT_LEN       <- 450L       # ~150 aa; full mcrA alpha is ~550 aa / ~1650 nt

N_THREADS        <- 8L

# Assembly-level filter (NULL = keep all)
ASSEMBLY_LEVELS  <- NULL       # e.g. c("Complete Genome","Chromosome")

# NCBI credentials
ENTREZ_EMAIL     <- "vasiliad@gmail.com"
rentrez::set_entrez_key("f7af40e4fb97a735bae78c1870be29914a08")

# Methanogenic archaeal lineages to query (NCBI taxonomy IDs)
METHANOGEN_TAXA <- list(
  Methanobacteria         = 183925L,
  Methanococci            = 183939L,
  Methanomicrobia         = 224756L,
  Methanopyri             = 183967L,
  Methanomassiliicoccales = 1458055L,
  Verstraetearchaeota     = 1986330L
)

NCBI_DATASETS_URL <- "https://api.ncbi.nlm.nih.gov/datasets/v2"

# mcrA annotation patterns (applied to CDS FASTA headers, case-insensitive)
# Matches:  [gene=mcrA], [gene=mcrA-I], [gene=mcrAI], etc.
#           methyl-coenzyme M reductase * alpha (various punctuation)
MCRA_HEADER_RE <- paste0(
  "(?i)(?:",
  "\\[gene=mcrA[^\\]]*\\]",
  "|methyl.{0,1}coenzyme.{0,3}M.{0,3}reductase[^\\n]{0,60}alpha",
  ")"
)

# ── 2.  Helper functions ───────────────────────────────────────────────────────

db_log <- function(msg, level = "INFO")
  cat(sprintf("%s  %-5s  %s\n", format(Sys.time(), "%H:%M:%S"), level, msg))

`%||%` <- function(x, y) if (!is.null(x) && length(x) > 0L) x else y

read_fasta <- function(path) {
  lines   <- readLines(path, warn = FALSE)
  hdr_idx <- grep("^>", lines)
  if (length(hdr_idx) == 0L) return(setNames(character(0), character(0)))
  starts  <- hdr_idx + 1L
  ends    <- c(hdr_idx[-1L] - 1L, length(lines))
  seqs    <- mapply(function(s, e) paste(lines[s:e], collapse = ""),
                    starts, ends, SIMPLIFY = TRUE)
  names(seqs) <- sub("^>\\s*", "", lines[hdr_idx])
  seqs
}

write_fasta <- function(seqs, path) {
  lines <- unlist(mapply(function(nm, sq) c(paste0(">", nm), sq),
                         names(seqs), seqs, SIMPLIFY = FALSE))
  writeLines(lines, path)
}

resolve_bin <- function(name) {
  exp <- path.expand(name)
  if (file.exists(exp) && file.access(exp, 1L) == 0L) return(exp)
  found <- Sys.which(name)
  if (nzchar(found)) return(found)
  shell_path <- tryCatch({
    raw <- system2("bash", c("-lc", shQuote(paste("which", name))),
                   stdout = TRUE, stderr = FALSE)
    if (length(raw) > 0L && nzchar(raw[1L])) raw[1L] else character(0)
  }, warning = function(w) character(0), error = function(e) character(0))
  if (length(shell_path) > 0L && file.exists(shell_path) &&
      file.access(shell_path, 1L) == 0L) return(shell_path)
  roots <- c("~/anaconda3/bin","~/miniconda3/bin","~/mambaforge/bin",
             "~/opt/anaconda3/bin","~/opt/miniconda3/bin",
             "/usr/local/bin","/opt/homebrew/bin",
             "~/ncbi","~/.local/bin","/usr/bin")
  for (r in roots) {
    p <- file.path(path.expand(r), name)
    if (file.exists(p) && file.access(p, 1L) == 0L) return(p)
  }
  name
}

run_bin <- function(bin, args, stdout = FALSE, stderr = FALSE) {
  tryCatch(
    system2(bin, args = args, stdout = stdout, stderr = stderr),
    warning = function(w) if (stdout || stderr) character(0) else 127L,
    error   = function(e) if (stdout || stderr) character(0) else 127L
  )
}

# ── enrich_lineages ────────────────────────────────────────────────────────────
enrich_lineages <- function(df) {
  LINEAGE_RANKS <- c("domain", "phylum", "class", "order", "family", "genus")
  if ("ncbi_domain" %in% names(df) &&
      sum(!is.na(df$ncbi_domain)) > nrow(df) * 0.30) {
    db_log("  Lineage columns already populated — skipping Entrez enrichment")
    return(df)
  }
  unique_taxids <- unique(na.omit(as.integer(df$tax_id)))
  if (length(unique_taxids) == 0L) {
    db_log("  No valid tax_ids — cannot enrich lineages", "WARN"); return(df)
  }
  db_log(sprintf("  Fetching lineages for %d unique tax IDs via Entrez XML …",
                 length(unique_taxids)))

  parse_taxonomy_xml <- function(xml_text, req_ids) {
    flat   <- gsub("[[:space:]]+", " ", xml_text)
    chunks <- strsplit(flat, "</LineageEx>", fixed = TRUE)[[1]]
    result <- list()
    for (chunk in chunks) {
      lin_pos <- regexpr("<LineageEx>", chunk, fixed = TRUE)
      if (lin_pos[1L] < 0L) next
      before_lin  <- substr(chunk, 1L, lin_pos[1L] - 1L)
      lin_content <- substr(chunk, lin_pos[1L] + nchar("<LineageEx>"), nchar(chunk))
      tids_before <- regmatches(before_lin,
        gregexpr("(?<=<TaxId>)\\d+(?=</TaxId>)", before_lin, perl = TRUE))[[1]]
      owner <- intersect(tids_before, req_ids)
      if (length(owner) == 0L) next
      owner <- owner[1L]
      sci <- regmatches(lin_content,
        gregexpr("(?<=<ScientificName>)[^<]+(?=</ScientificName>)",
                 lin_content, perl = TRUE))[[1]]
      rnk <- regmatches(lin_content,
        gregexpr("(?<=<Rank>)[^<]+(?=</Rank>)", lin_content, perl = TRUE))[[1]]
      lin <- list()
      for (i in seq_len(min(length(sci), length(rnk)))) {
        rk <- tolower(trimws(rnk[i]))
        if (rk %in% LINEAGE_RANKS) lin[[rk]] <- trimws(sci[i])
      }
      result[[owner]] <- lin
    }
    result
  }

  lineage_map <- list()
  batches     <- split(unique_taxids, ceiling(seq_along(unique_taxids) / 200L))
  for (b in seq_along(batches)) {
    ids      <- as.character(batches[[b]])
    xml_text <- tryCatch(
      rentrez::entrez_fetch(db = "taxonomy", id = ids, rettype = "xml"),
      error = function(e) {
        db_log(sprintf("  Lineage batch %d/%d failed: %s",
                       b, length(batches), e$message), "WARN"); NULL
      }
    )
    Sys.sleep(0.12)
    if (is.null(xml_text) || !nzchar(xml_text)) next
    lineage_map <- c(lineage_map, parse_taxonomy_xml(xml_text, ids))
    if (b %% 5L == 0L || b == length(batches))
      db_log(sprintf("  Lineage fetch: %d / %d batches (%d mapped)",
                     b, length(batches), length(lineage_map)))
  }

  for (rk in LINEAGE_RANKS) {
    col_nm <- paste0("ncbi_", rk)
    df[[col_nm]] <- vapply(as.character(df$tax_id), function(tid) {
      if (is.na(tid) || !tid %in% names(lineage_map)) return(NA_character_)
      lineage_map[[tid]][[rk]] %||% NA_character_
    }, character(1L))
  }
  db_log(sprintf("  Lineage enrichment done: %d / %d genomes have phylum",
                 sum(!is.na(df$ncbi_phylum)), nrow(df)))
  df
}

# ── fetch_ncbi_genome_reports ──────────────────────────────────────────────────
fetch_ncbi_genome_reports <- function(taxon_id, taxon_name, page_size = 1000L) {
  url        <- sprintf("%s/genome/taxon/%d/dataset_report", NCBI_DATASETS_URL, taxon_id)
  all_rows   <- list()
  page_token <- NULL
  repeat {
    params <- list("filters.reference_only" = "false",
                   "filters.exclude_paired_reports" = "false",
                   "page_size" = as.character(page_size))
    if (!is.null(page_token)) params[["page_token"]] <- page_token
    resp <- tryCatch(
      httr::GET(url, query = params,
                httr::add_headers("Accept" = "application/json"),
                httr::timeout(60)),
      error = function(e) NULL
    )
    if (is.null(resp) || httr::status_code(resp) != 200L) {
      db_log(sprintf("  [%s] API request failed (status %s)", taxon_name,
                     if (is.null(resp)) "TIMEOUT" else httr::status_code(resp)), "WARN")
      break
    }
    body    <- jsonlite::fromJSON(httr::content(resp, as = "text", encoding = "UTF-8"),
                                  simplifyVector = FALSE)
    reports <- body[["reports"]]
    if (length(reports) == 0L) break
    for (rpt in reports) {
      asm <- rpt[["assembly_info"]]; org <- rpt[["organism"]]
      lineage <- list()
      cls <- rpt[["taxonomy"]][["classification"]]
      if (!is.null(cls))
        for (node in cls) lineage[[tolower(node[["rank"]])]] <- node[["name"]]
      all_rows <- c(all_rows, list(data.frame(
        accession      = rpt[["accession"]]              %||% NA_character_,
        organism_name  = org[["organism_name"]]          %||% NA_character_,
        tax_id         = org[["tax_id"]]                 %||% NA_integer_,
        assembly_level = asm[["assembly_level"]]         %||% NA_character_,
        assembly_name  = asm[["assembly_name"]]          %||% NA_character_,
        ncbi_domain    = lineage[["domain"]]             %||% NA_character_,
        ncbi_phylum    = lineage[["phylum"]]             %||% NA_character_,
        ncbi_class     = lineage[["class"]]              %||% NA_character_,
        ncbi_order     = lineage[["order"]]              %||% NA_character_,
        ncbi_family    = lineage[["family"]]             %||% NA_character_,
        ncbi_genus     = lineage[["genus"]]              %||% NA_character_,
        ncbi_species   = org[["organism_name"]]          %||% NA_character_,
        query_group    = taxon_name,
        stringsAsFactors = FALSE
      )))
    }
    page_token <- body[["next_page_token"]]
    if (is.null(page_token) || !nzchar(page_token)) break
    Sys.sleep(0.2)
  }
  if (length(all_rows) == 0L) return(data.frame())
  do.call(rbind, all_rows)
}

# ── CDS download helpers (datasets CLI + FTP fallback) ────────────────────────
load_assembly_ftp_map <- function() {
  if (file.exists(ASM_FTP_CACHE) &&
      as.numeric(difftime(Sys.time(), file.mtime(ASM_FTP_CACHE), units = "days")) < 7) {
    ftp_map <- readRDS(ASM_FTP_CACHE)
    # Validate: cache must contain GCF_ entries (RefSeq); if not, it was built
    # before the dual GenBank+RefSeq download was added and must be rebuilt.
    if (!any(startsWith(names(ftp_map), "GCF_"))) {
      db_log("  FTP map cache has no GCF_ entries — rebuilding …", "WARN")
      unlink(ASM_FTP_CACHE)
    } else {
      db_log(sprintf("  Assembly FTP map loaded from cache (%d entries)", length(ftp_map)))
      return(ftp_map)
    }
  }
  asm_urls <- c(
    GenBank = "https://ftp.ncbi.nlm.nih.gov/genomes/genbank/archaea/assembly_summary.txt",
    RefSeq  = "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/archaea/assembly_summary.txt"
  )
  parse_asm_summary <- function(url, label) {
    db_log(sprintf("  Downloading %s archaea assembly_summary.txt …", label))
    tmp <- tempfile(fileext = ".txt")
    res <- tryCatch(httr::GET(url, httr::write_disk(tmp, overwrite = TRUE),
                               httr::timeout(300)), error = function(e) NULL)
    if (is.null(res) || httr::status_code(res) != 200L) {
      db_log(sprintf("  %s download failed", label), "WARN"); unlink(tmp); return(NULL)
    }
    lines   <- readLines(tmp, warn = FALSE); unlink(tmp)
    hdr_idx <- grep("assembly_accession", lines, fixed = TRUE)[1L]
    if (is.na(hdr_idx)) return(NULL)
    lines[hdr_idx] <- sub("^#\\s*", "", lines[hdr_idx])
    df <- tryCatch(
      read.table(text = paste(lines[hdr_idx:length(lines)], collapse = "\n"),
                 sep = "\t", header = TRUE, quote = "", comment.char = "",
                 check.names = FALSE, stringsAsFactors = FALSE, fill = TRUE),
      error = function(e) NULL
    )
    if (is.null(df) || !all(c("assembly_accession","ftp_path") %in% names(df))) return(NULL)
    df
  }
  dfs <- Filter(Negate(is.null),
                lapply(seq_along(asm_urls), function(i)
                  parse_asm_summary(asm_urls[i], names(asm_urls)[i])))
  if (length(dfs) == 0L) return(setNames(character(0), character(0)))
  df_all  <- do.call(rbind, dfs)
  valid   <- df_all$ftp_path != "na" & nzchar(df_all$ftp_path)
  acc_vec <- df_all$assembly_accession[valid]
  ftp_vec <- sub("/$", "", sub("^ftp://", "https://", df_all$ftp_path[valid]))
  dup     <- duplicated(acc_vec)
  ftp_map <- setNames(ftp_vec[!dup], acc_vec[!dup])
  saveRDS(ftp_map, ASM_FTP_CACHE)
  db_log(sprintf("  Assembly FTP map built: %d entries", length(ftp_map)))
  ftp_map
}

datasets_bin <- resolve_bin("datasets")
datasets_ok  <- FALSE   # set below after binary check

download_cds_datasets <- function(accession) {
  out_fa  <- file.path(GENOME_CDS_DIR, paste0(accession, ".fna"))
  out_zip <- file.path(GENOME_CDS_DIR, paste0(accession, ".zip"))
  if (file.exists(out_fa) && file.info(out_fa)$size > 0L) return(out_fa)
  ret <- run_bin(datasets_bin,
                 args = c("download","genome","accession", accession,
                           "--include","cds","--filename", out_zip,
                           "--no-progressbar"))
  if (ret != 0L || !file.exists(out_zip)) return(NULL)
  tmp_dir   <- file.path(GENOME_CDS_DIR, paste0("tmp_", accession))
  unzip(out_zip, exdir = tmp_dir)
  cds_files <- list.files(tmp_dir, pattern = "cds_from_genomic\\.fna(\\.gz)?$",
                           recursive = TRUE, full.names = TRUE)
  if (length(cds_files) == 0L) { unlink(tmp_dir, recursive = TRUE); return(NULL) }
  src <- cds_files[1L]
  if (grepl("\\.gz$", src)) {
    con_in  <- gzcon(file(src, "rb")); con_out <- file(out_fa, "wb")
    writeBin(readBin(con_in, "raw", n = 200e6L), con_out)
    close(con_in); close(con_out)
  } else file.copy(src, out_fa)
  unlink(tmp_dir, recursive = TRUE)
  if (file.exists(out_fa) && file.info(out_fa)$size > 0L) out_fa else NULL
}

download_cds_ftp <- function(accession, ftp_map = NULL) {
  out_fa <- file.path(GENOME_CDS_DIR, paste0(accession, ".fna"))
  if (file.exists(out_fa) && file.info(out_fa)$size > 0L) return(out_fa)
  ftp_base <- NA_character_
  if (!is.null(ftp_map) && accession %in% names(ftp_map))
    ftp_base <- ftp_map[[accession]]
  if (is.na(ftp_base)) {
    asm_name <- tryCatch({
      base_acc <- sub("\\.\\d+$", "", accession)   # strip version (e.g. .1)
      s <- rentrez::entrez_search(db = "assembly", term = accession,
                                  use_history = FALSE)
      Sys.sleep(0.11)
      if (length(s$ids) == 0L) {
        s <- rentrez::entrez_search(db = "assembly", term = base_acc,
                                    use_history = FALSE)
        Sys.sleep(0.11)
      }
      if (length(s$ids) == 0L) stop("no hits")
      summ  <- rentrez::entrez_summary(db = "assembly", id = s$ids[1L])
      Sys.sleep(0.11)
      aname <- summ[["assemblyname"]]
      if (!is.null(aname) && nzchar(aname)) aname else stop("empty name")
    }, error = function(e) NA_character_)
    if (!is.na(asm_name)) {
      base <- sub("\\..*$", "", accession)
      pts  <- strsplit(base, "_", fixed = TRUE)[[1L]]
      num  <- sprintf("%09s", pts[2L])
      ftp_base <- sprintf(
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/%s/%s/%s/%s/%s_%s",
        pts[1L], substr(num,1,3), substr(num,4,6), substr(num,7,9),
        accession, asm_name)
    }
  }
  if (is.na(ftp_base)) {
    db_log(sprintf("  [FTP] %s: FTP path unresolvable — skipping", accession), "WARN")
    return(NULL)
  }
  ftp_base  <- sub("/$", "", ftp_base)
  dir_name  <- basename(ftp_base)
  cds_url   <- paste0(ftp_base, "/", dir_name, "_cds_from_genomic.fna.gz")
  tmp_gz    <- tempfile(fileext = ".fna.gz")
  res <- tryCatch(httr::GET(cds_url, httr::write_disk(tmp_gz, overwrite = TRUE),
                             httr::timeout(180)), error = function(e) NULL)
  if (is.null(res) || httr::status_code(res) != 200L) {
    unlink(tmp_gz); return(NULL)
  }
  tryCatch({
    con_in  <- gzcon(file(tmp_gz, "rb")); con_out <- file(out_fa, "wb")
    writeBin(readBin(con_in, "raw", n = 200e6L), con_out)
    close(con_in); close(con_out)
  }, error = function(e) NULL)
  unlink(tmp_gz)
  if (file.exists(out_fa) && file.info(out_fa)$size > 0L) out_fa else NULL
}

# =============================================================================
# STEP 1 — Genome accessions + lineage enrichment
# =============================================================================
db_log("=== STEP 1: Fetching methanogen genome accessions ===")

for (d in c(CACHE_DIR, GENOME_CDS_DIR, MCRA_RAW_DIR))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

if (file.exists(ACCESSION_CACHE)) {
  db_log("  Loading accession list from cache …")
  genome_df <- readRDS(ACCESSION_CACHE)
  if (!"assembly_name" %in% names(genome_df)) {
    db_log("  Cache missing assembly_name — re-querying API …", "WARN")
    unlink(ACCESSION_CACHE)
    genome_df <- NULL
  }
}
if (!file.exists(ACCESSION_CACHE)) {
  db_log("  Querying NCBI Datasets API …")
  parts     <- mapply(fetch_ncbi_genome_reports,
                      taxon_id   = unlist(METHANOGEN_TAXA),
                      taxon_name = names(METHANOGEN_TAXA),
                      SIMPLIFY   = FALSE)
  genome_df <- dplyr::bind_rows(parts)
  genome_df <- genome_df[!duplicated(genome_df$accession), ]
  saveRDS(genome_df, ACCESSION_CACHE)
}

db_log(sprintf("  Total genomes: %d", nrow(genome_df)))
print(table(genome_df$assembly_level))

if (!is.null(ASSEMBLY_LEVELS)) {
  genome_df <- genome_df[genome_df$assembly_level %in% ASSEMBLY_LEVELS, ]
  db_log(sprintf("  After level filter: %d genomes", nrow(genome_df)))
}
if (nrow(genome_df) == 0L) stop("No genomes passed filters. Check ASSEMBLY_LEVELS.")

genome_df <- enrich_lineages(genome_df)
saveRDS(genome_df, ACCESSION_CACHE)

# =============================================================================
# STEP 2 — Download CDS sequences per genome
# =============================================================================
db_log("=== STEP 2: Downloading CDS sequences ===")

datasets_ok <- (run_bin(datasets_bin, "--version") == 0L)
db_log(sprintf("  datasets CLI usable: %s", datasets_ok))
if (!datasets_ok) {
  db_log("  → falling back to NCBI FTP direct download", "WARN")
  ftp_map <- load_assembly_ftp_map()
} else {
  ftp_map <- NULL
}

genome_df$cds_path <- NA_character_
for (i in seq_len(nrow(genome_df))) {
  acc <- genome_df$accession[i]
  p   <- if (datasets_ok) download_cds_datasets(acc) else download_cds_ftp(acc, ftp_map)
  genome_df$cds_path[i] <- if (!is.null(p)) p else NA_character_
  if (i %% 100L == 0L)
    db_log(sprintf("  [%d / %d]  downloaded=%d  failed=%d",
                   i, nrow(genome_df),
                   sum(!is.na(genome_df$cds_path[seq_len(i)])),
                   sum( is.na(genome_df$cds_path[seq_len(i)]))))
}

genome_df_ok <- genome_df[!is.na(genome_df$cds_path), ]
db_log(sprintf("  CDS downloaded: %d / %d genomes",
               nrow(genome_df_ok), nrow(genome_df)))

# =============================================================================
# STEP 3 — Extract mcrA-annotated CDS by header pattern
# =============================================================================
# NCBI CDS FASTA headers contain annotation tags such as:
#   [gene=mcrA] [protein=methyl-coenzyme M reductase subunit alpha]
# We extract all CDS matching MCRA_HEADER_RE and write one FASTA per genome.
# =============================================================================
db_log("=== STEP 3: Extracting mcrA-annotated CDS by header ===")

extract_mcra_cds <- function(acc, cds_fa) {
  out_fa <- file.path(MCRA_RAW_DIR, paste0(acc, "_mcra_raw.fna"))
  if (file.exists(out_fa) && file.info(out_fa)$size > 0L) return(out_fa)

  seqs <- tryCatch(read_fasta(cds_fa), error = function(e) NULL)
  if (is.null(seqs) || length(seqs) == 0L) return(NULL)

  # Match header against mcrA annotation patterns
  is_mcra  <- grepl(MCRA_HEADER_RE, names(seqs), perl = TRUE)
  mcra_seqs <- seqs[is_mcra]
  if (length(mcra_seqs) == 0L) return(NULL)

  # Rename: prepend accession so every sequence has a globally unique ID
  names(mcra_seqs) <- paste0(acc, "|", names(mcra_seqs))
  write_fasta(mcra_seqs, out_fa)
  out_fa
}

genome_df_ok$mcra_raw_path <- NA_character_
n_total_mcra <- 0L
for (i in seq_len(nrow(genome_df_ok))) {
  acc <- genome_df_ok$accession[i]
  p   <- tryCatch(extract_mcra_cds(acc, genome_df_ok$cds_path[i]),
                  error = function(e) NULL)
  genome_df_ok$mcra_raw_path[i] <- if (!is.null(p)) p else NA_character_
}

mcra_files <- genome_df_ok$mcra_raw_path[!is.na(genome_df_ok$mcra_raw_path)]
n_genomes_with_mcra <- length(mcra_files)
all_mcra_raw <- unlist(lapply(mcra_files, read_fasta), use.names = TRUE)
db_log(sprintf("  Genomes with ≥1 mcrA CDS: %d / %d",
               n_genomes_with_mcra, nrow(genome_df_ok)))
db_log(sprintf("  Total mcrA CDS extracted : %d", length(all_mcra_raw)))

if (length(all_mcra_raw) == 0L)
  stop("No mcrA-annotated CDS found. Check MCRA_HEADER_RE and CDS downloads.")

# =============================================================================
# STEP 4 — HMMER completeness filter
# =============================================================================
# Translate each extracted CDS to protein (already in frame 1) and run
# hmmsearch against TIGR03256.hmm.  Keep sequences where the HMM domain
# alignment covers >= HMM_MIN_COVERAGE of the profile length.
# Trim the corresponding nucleotide sequence to the envelope coordinates.
# =============================================================================
db_log("=== STEP 4: HMMER completeness filter (TIGR03256) ===")

if (!file.exists(HMM_PROFILE))
  stop("Missing HMM profile: ", HMM_PROFILE,
       "\nDownload from https://tigrfams.jcvi.org/cgi-bin/HmmReportPage.cgi?acc=TIGR03256")

hmmsearch_bin <- resolve_bin("hmmsearch")
hmmpress_bin  <- resolve_bin("hmmpress")

if (run_bin(hmmsearch_bin, "--version") == 127L)
  stop("hmmsearch not found. Install HMMER3: conda install -c bioconda hmmer")

# Press the HMM if needed
if (!file.exists(paste0(HMM_PROFILE, ".h3i"))) {
  db_log("  Running hmmpress …")
  run_bin(hmmpress_bin, args = HMM_PROFILE)
}

# Get HMM profile length from the profile header
hmm_hdr <- readLines(HMM_PROFILE, n = 30L)
hmm_len  <- as.integer(
  sub("LENG\\s+", "", hmm_hdr[grep("^LENG", hmm_hdr)][1L]))
db_log(sprintf("  TIGR03256 profile length: %d aa", hmm_len))

# ── Translate CDS sequences ────────────────────────────────────────────────────
# NCBI CDS sequences are in reading frame +1 and include the stop codon.
# We translate directly, replacing stop codons (*) with X.
translate_cds <- function(nt_seq) {
  # Trim to a multiple of 3 codons (drop partial codon at end)
  trim_len <- floor(nchar(nt_seq) / 3L) * 3L
  nt_trim  <- substr(nt_seq, 1L, trim_len)
  if (nchar(nt_trim) < 3L) return(NA_character_)
  aa <- tryCatch(
    suppressWarnings(
      as.character(Biostrings::translate(
        Biostrings::DNAString(nt_trim), if.fuzzy.codon = "X"))
    ),
    error = function(e) NA_character_
  )
  if (is.na(aa)) return(NA_character_)
  # Replace stop codon symbol and trim terminal stop
  aa <- gsub("\\*", "X", aa)
  aa <- sub("X$", "", aa)         # remove trailing stop (common in NCBI CDS)
  aa
}

db_log(sprintf("  Translating %d CDS sequences …", length(all_mcra_raw)))
aa_seqs <- vapply(all_mcra_raw, translate_cds, character(1L))
valid_aa <- !is.na(aa_seqs) & nchar(aa_seqs) >= 50L
aa_seqs  <- aa_seqs[valid_aa]
db_log(sprintf("  Translatable sequences: %d", length(aa_seqs)))

# ── Write protein FASTA and run hmmsearch ─────────────────────────────────────
# hmmsearch truncates FASTA headers at the first whitespace, so the seq_id in
# the domain table would not match the full original names (which contain
# annotation tags with spaces: "[gene=mcrA] [protein=...]").
# Solution: assign short numeric IDs for the protein FASTA, then map back
# after parsing the domain table using safe_to_orig.
safe_ids    <- sprintf("s%08d", seq_along(aa_seqs))
safe_to_orig <- setNames(names(aa_seqs), safe_ids)   # safe_id → original name
aa_seqs_safe <- setNames(aa_seqs, safe_ids)

tmp_prot_fa <- tempfile(fileext = ".faa")
tmp_dom_out <- tempfile(fileext = ".domtbl")
write_fasta(aa_seqs_safe, tmp_prot_fa)

db_log("  Running hmmsearch …")
run_bin(hmmsearch_bin,
        args = c("--domtblout", tmp_dom_out,
                 "--noali", "-E", "1e-5",
                 "--cpu", as.character(N_THREADS),
                 HMM_PROFILE, tmp_prot_fa))

# ── Parse domain table ─────────────────────────────────────────────────────────
dom_lines <- readLines(tmp_dom_out, warn = FALSE)
dom_lines <- dom_lines[!startsWith(dom_lines, "#") & nzchar(trimws(dom_lines))]

dom_df <- do.call(rbind, lapply(dom_lines, function(l) {
  f <- strsplit(trimws(l), "\\s+")[[1L]]
  if (length(f) < 22L) return(NULL)
  data.frame(
    seq_id   = f[1L],
    hmm_from = as.integer(f[16L]),
    hmm_to   = as.integer(f[17L]),
    env_from = as.integer(f[20L]),   # AA envelope start (1-based)
    env_to   = as.integer(f[21L]),   # AA envelope end   (1-based)
    score    = as.numeric(f[14L]),
    stringsAsFactors = FALSE
  )
}))

unlink(c(tmp_prot_fa, tmp_dom_out))

# Translate safe IDs back to original sequence names
if (!is.null(dom_df) && nrow(dom_df) > 0L)
  dom_df$seq_id <- safe_to_orig[dom_df$seq_id]

if (is.null(dom_df) || nrow(dom_df) == 0L)
  stop("hmmsearch returned no hits. Check HMM_PROFILE and input sequences.")

# Keep best-scoring domain hit per sequence; compute HMM coverage
dom_df$hmm_cov <- (dom_df$hmm_to - dom_df$hmm_from + 1L) / hmm_len
dom_df         <- dom_df[order(dom_df$seq_id, -dom_df$score), ]
dom_df         <- dom_df[!duplicated(dom_df$seq_id), ]

db_log(sprintf("  hmmsearch hits: %d", nrow(dom_df)))
db_log(sprintf("  Coverage distribution (quantiles):"))
print(round(quantile(dom_df$hmm_cov, c(0.1, 0.25, 0.5, 0.75, 0.9)), 3))

pass_hmm <- dom_df[dom_df$hmm_cov >= HMM_MIN_COVERAGE, ]
db_log(sprintf("  Sequences passing >= %.0f%% HMM coverage: %d / %d",
               HMM_MIN_COVERAGE * 100, nrow(pass_hmm), nrow(dom_df)))

if (nrow(pass_hmm) == 0L)
  stop("No sequences passed the HMM completeness filter. Lower HMM_MIN_COVERAGE?")

# =============================================================================
# STEP 5 — Trim nucleotide sequences to HMM envelope coordinates
# =============================================================================
# The hmmsearch domain table reports envelope coordinates in AA space.
# Convert to nucleotide coordinates and trim the original CDS.
# =============================================================================
db_log("=== STEP 5: Trimming sequences to HMM envelope ===")

trim_to_envelope <- function(seq_id, env_from_aa, env_to_aa) {
  nt_seq <- all_mcra_raw[[seq_id]]
  if (is.null(nt_seq) || !nzchar(nt_seq)) return(NULL)

  # AA coordinates → nucleotide coordinates (frame +1, 1-based)
  nt_start <- (env_from_aa - 1L) * 3L + 1L
  nt_end   <- env_to_aa * 3L
  nt_end   <- min(nt_end, nchar(nt_seq))

  trimmed <- substr(nt_seq, nt_start, nt_end)
  if (nchar(trimmed) < MIN_NT_LEN) return(NULL)
  trimmed
}

trimmed_seqs <- list()
for (k in seq_len(nrow(pass_hmm))) {
  sid <- pass_hmm$seq_id[k]
  tr  <- trim_to_envelope(sid, pass_hmm$env_from[k], pass_hmm$env_to[k])
  if (!is.null(tr)) trimmed_seqs[[sid]] <- tr
}
trimmed_seqs <- unlist(trimmed_seqs)

db_log(sprintf("  Sequences after envelope trim: %d", length(trimmed_seqs)))

# =============================================================================
# STEP 6 — Deduplicate by exact nucleotide sequence
# =============================================================================
db_log("=== STEP 6: Deduplicating by exact sequence ===")

n_before   <- length(trimmed_seqs)
keep_mask  <- !duplicated(toupper(trimmed_seqs))
uniq_seqs  <- trimmed_seqs[keep_mask]
db_log(sprintf("  %d raw → %d unique sequences", n_before, length(uniq_seqs)))

# =============================================================================
# STEP 7 — Build DADA2 taxonomy headers
# =============================================================================
# Header format: >Domain;Phylum;Class;Order;Family;Genus;Species;Accession
# The genome accession (GCA_/GCF_) is appended as the final field so every
# database entry can be traced back to its source genome.
# Missing ranks propagate the last known rank with "_unclassified" appended.
# =============================================================================
db_log("=== STEP 7: Building DADA2 taxonomy headers ===")

RANKS <- c("domain","phylum","class","order","family","genus","species")

# Build accession → taxonomy lookup from genome_df
tax_lookup <- setNames(
  lapply(seq_len(nrow(genome_df)), function(i) {
    r <- genome_df[i, ]
    list(
      domain  = r$ncbi_domain  %||% NA_character_,
      phylum  = r$ncbi_phylum  %||% NA_character_,
      class   = r$ncbi_class        %||% NA_character_,
      order   = r$ncbi_order        %||% NA_character_,
      family  = r$ncbi_family       %||% NA_character_,
      genus   = r$ncbi_genus        %||% NA_character_,
      species = r$ncbi_species      %||% NA_character_
    )
  }),
  genome_df$accession
)

# For missing ranks, propagate parent rank with "_unclassified" suffix
fill_lineage <- function(lin_list) {
  vals <- character(length(RANKS))
  names(vals) <- RANKS
  last_known <- NA_character_
  for (rk in RANKS) {
    v <- lin_list[[rk]]
    if (!is.na(v) && nzchar(trimws(v))) {
      vals[rk] <- trimws(v)
      last_known <- trimws(v)
    } else {
      vals[rk] <- if (!is.na(last_known)) paste0(last_known, "_unclassified")
                  else "Unclassified"
    }
  }
  vals
}

dada2_seqs <- character(length(uniq_seqs))
for (i in seq_along(uniq_seqs)) {
  nm  <- names(uniq_seqs)[i]
  acc <- sub("\\|.*$", "", nm)                   # genome accession

  lin <- if (acc %in% names(tax_lookup)) {
    fill_lineage(tax_lookup[[acc]])
  } else {
    setNames(rep("Unclassified", length(RANKS)), RANKS)
  }

  tax_hdr <- paste(c(lin, acc), collapse = ";")
  dada2_seqs[i] <- uniq_seqs[i]
  names(dada2_seqs)[i] <- tax_hdr
}

# =============================================================================
# STEP 8 — Write output
# =============================================================================
db_log("=== STEP 8: Writing output ===")

write_fasta(dada2_seqs, OUT_FASTA)
db_log(sprintf("  Database written: %s  (%d sequences)", OUT_FASTA, length(dada2_seqs)))

# Per-rank sequence count statistics (ranks 1–7; position 8 is genome accession)
stats_rows <- lapply(seq_along(RANKS), function(ri) {
  rk <- RANKS[ri]
  tax_vals <- sapply(strsplit(names(dada2_seqs), ";"), `[`, ri)
  n_unclass <- sum(grepl("_unclassified|^Unclassified$", tax_vals))
  data.frame(
    rank          = rk,
    n_sequences   = length(dada2_seqs),
    n_classified  = length(dada2_seqs) - n_unclass,
    n_unclassified= n_unclass,
    n_unique_taxa = length(unique(tax_vals[!grepl("_unclassified|^Unclassified$", tax_vals)])),
    stringsAsFactors = FALSE
  )
})
stats_df <- do.call(rbind, stats_rows)
write.table(stats_df, OUT_STATS, sep = "\t", row.names = FALSE, quote = FALSE)

# Console summary
cat("\n")
cat(strrep("═", 65), "\n")
cat("  mcrA NCBI-ANNOTATION DATABASE — SUMMARY\n")
cat(strrep("═", 65), "\n")
cat(sprintf("  %-14s %8s %12s %12s\n",
            "Rank", "Seqs", "Classified", "Unique taxa"))
cat(strrep("─", 65), "\n")
for (i in seq_len(nrow(stats_df))) {
  r <- stats_df[i, ]
  cat(sprintf("  %-14s %8d %12d %12d\n",
              r$rank, r$n_sequences, r$n_classified, r$n_unique_taxa))
}
cat(strrep("═", 65), "\n")
cat(sprintf("\n  Output : %s\n", OUT_FASTA))
cat(sprintf("  Stats  : %s\n\n", OUT_STATS))

db_log("=== Done ===")

# =============================================================================
# STEP 9 — Build yang_2014_taxize.fasta
# =============================================================================
# Joins mcrAtemplate.fasta sequences with NCBI-fetched lineages to produce a
# DADA2-format Yang 2014 database with complete, up-to-date taxonomy.
# Results cached in yang_lineage_cache.rds; delete to force re-fetch.
# =============================================================================
YANG_TEMPLATE_FA  <- "mcrAtemplate.fasta"
YANG_TAXONOMY_TXT <- "Yang_etal_2014_tax4mcrA.taxonomy.txt"
YANG_TAXIZE_FA    <- "yang_2014_taxize.fasta"
YANG_CACHE_RDS    <- "yang_lineage_cache.rds"

if (file.exists(YANG_TEMPLATE_FA) && file.exists(YANG_TAXONOMY_TXT)) {
  db_log("=== STEP 9: Building yang_2014_taxize.fasta ===")

  # ── parse_taxonomy_xml_local helper ────────────────────────────────────────
  .parse_yang_xml <- function(xml_text, req_ids) {
    LRANKS <- c("superkingdom","phylum","class","order","family","genus","species")
    flat   <- gsub("[[:space:]]+", " ", xml_text)
    chunks <- strsplit(flat, "</LineageEx>", fixed = TRUE)[[1]]
    result <- list()
    for (chunk in chunks) {
      lin_pos <- regexpr("<LineageEx>", chunk, fixed = TRUE)
      if (lin_pos[1L] < 0L) next
      before_lin  <- substr(chunk, 1L, lin_pos[1L] - 1L)
      lin_content <- substr(chunk, lin_pos[1L] + nchar("<LineageEx>"), nchar(chunk))
      tids_before <- regmatches(before_lin,
        gregexpr("(?<=<TaxId>)\\d+(?=</TaxId>)", before_lin, perl = TRUE))[[1]]
      owner <- intersect(tids_before, req_ids)
      if (length(owner) == 0L) next
      owner <- owner[1L]
      sci <- regmatches(lin_content,
        gregexpr("(?<=<ScientificName>)[^<]+(?=</ScientificName>)",
                 lin_content, perl = TRUE))[[1]]
      rnk <- regmatches(lin_content,
        gregexpr("(?<=<Rank>)[^<]+(?=</Rank>)", lin_content, perl = TRUE))[[1]]
      lin <- list()
      for (i in seq_len(min(length(sci), length(rnk)))) {
        rk <- tolower(trimws(rnk[i]))
        if (rk == "domain") rk <- "superkingdom"
        if (rk %in% LRANKS) lin[[rk]] <- trimws(sci[i])
      }
      leaf_sci <- regmatches(xml_text,
        regexpr(sprintf("(?<=<TaxId>%s</TaxId>[^<]{0,200}<ScientificName>)[^<]+",
                        owner), xml_text, perl = TRUE))
      if (length(leaf_sci) > 0L && !("species" %in% names(lin)))
        lin[["species"]] <- trimws(leaf_sci[1L])
      result[[owner]] <- lin
    }
    result
  }

  # ── Build cache if absent ───────────────────────────────────────────────────
  if (!file.exists(YANG_CACHE_RDS)) {
    db_log("  Fetching NCBI taxonomy lineages for Yang 2014 …")

    yang_raw_df <- read.table(YANG_TAXONOMY_TXT, sep = "\t",
                               col.names = c("accession","tax_string"),
                               stringsAsFactors = FALSE) %>%
      transform(tax_string = sub(";$", "", tax_string))

    yang_raw_df$best_name <- sapply(
      strsplit(yang_raw_df$tax_string, ";", fixed = TRUE),
      function(parts) {
        for (p in rev(parts)) {
          p <- trimws(gsub("_", " ", p))
          if (nzchar(p)) return(p)
        }
        NA_character_
      }
    )

    unique_names <- sort(unique(na.omit(yang_raw_df$best_name)))
    db_log(sprintf("  Unique organism names: %d", length(unique_names)))

    name_to_taxid <- setNames(character(length(unique_names)), unique_names)
    for (i in seq_along(unique_names)) {
      nm <- unique_names[i]
      res <- tryCatch(
        rentrez::entrez_search(db = "taxonomy",
                               term = sprintf('"%s"[Scientific Name]', nm),
                               retmax = 1L),
        error = function(e) NULL
      )
      Sys.sleep(0.12)
      if ((is.null(res) || length(res$ids) == 0L) &&
          !grepl("^uncultured|^unclassified|^environmental", nm, ignore.case = TRUE)) {
        res <- tryCatch(
          rentrez::entrez_search(db = "taxonomy",
                                 term = sprintf('%s[All Fields]', nm),
                                 retmax = 1L),
          error = function(e) NULL
        )
        Sys.sleep(0.12)
      }
      if (!is.null(res) && length(res$ids) > 0L)
        name_to_taxid[nm] <- res$ids[1L]
      if (i %% 25L == 0L || i == length(unique_names))
        db_log(sprintf("  Name search %d/%d (resolved: %d)",
                       i, length(unique_names), sum(nzchar(name_to_taxid))))
    }

    unique_taxids <- unique(name_to_taxid[nzchar(name_to_taxid)])
    taxid_lineage <- list()
    for (batch in split(unique_taxids, ceiling(seq_along(unique_taxids) / 200L))) {
      xml_text <- tryCatch(
        rentrez::entrez_fetch(db = "taxonomy", id = as.character(batch),
                               rettype = "xml"),
        error = function(e) NULL
      )
      Sys.sleep(0.12)
      if (!is.null(xml_text) && nzchar(xml_text))
        taxid_lineage <- c(taxid_lineage,
                           .parse_yang_xml(xml_text, as.character(batch)))
    }
    db_log(sprintf("  Lineages fetched: %d", length(taxid_lineage)))

    acc_to_taxid <- setNames(name_to_taxid[yang_raw_df$best_name],
                             yang_raw_df$accession)
    saveRDS(list(acc_to_taxid  = acc_to_taxid,
                 taxid_lineage = taxid_lineage,
                 name_to_taxid = name_to_taxid),
            YANG_CACHE_RDS)
    db_log(sprintf("  Cache saved → %s", YANG_CACHE_RDS))
  } else {
    db_log(sprintf("  Loading Yang cache from %s …", YANG_CACHE_RDS))
  }

  yang_c        <- readRDS(YANG_CACHE_RDS)
  acc_to_taxid  <- yang_c$acc_to_taxid
  taxid_lineage <- yang_c$taxid_lineage

  # ── Build yang_2014_taxize.fasta ────────────────────────────────────────────
  if (!file.exists(YANG_TAXIZE_FA)) {
    LRANKS_MAP <- c(superkingdom="Domain", phylum="Phylum", class="Class",
                    order="Order", family="Family", genus="Genus", species="Species")

    yang_raw_lines <- readLines(YANG_TAXONOMY_TXT)
    yang_map <- setNames(
      trimws(sub(";$", "", sub(".*\t", "", yang_raw_lines))),
      trimws(sub("\t.*", "", yang_raw_lines))
    )

    pad_to_7 <- function(s) {
      p <- strsplit(s, ";", fixed = TRUE)[[1L]]
      n <- length(p)
      if (n >= 7L) return(paste(p[1:7], collapse = ";"))
      if (n == 0L) return(paste(rep("Unclassified", 7L), collapse = ";"))
      paste(c(p[1L], rep("Unclassified", 7L - n), p[2L:n]), collapse = ";")
    }

    # ── Rank-specific nomenclature normalisation ───────────────────────────────
    # Each rename is applied only to the rank column where it belongs, preventing
    # cross-rank contamination (e.g. "Euryarchaeota" at Class due to depth shifts).
    # Family rename precedes Genus rename to avoid the substring collision:
    # "Methanosaeta" is a prefix of "Methanosaetaceae", so replacing Genus/Species
    # first would corrupt Family "Methanosaetaceae" → "Methanothrixceae".
    normalise_yang_tax <- function(tax_str) {
      p <- strsplit(tax_str, ";", fixed = TRUE)[[1L]]
      if (length(p) < 7L) p <- c(p, rep("Unclassified", 7L - length(p)))
      # All ranks: underscore → space, trim whitespace
      p <- trimws(gsub("_", " ", p))
      # Phylum (2): Euryarchaeota → Methanobacteriota (GTDB reclassification)
      p[2L] <- gsub("Euryarchaeota",    "Methanobacteriota", p[2L], fixed = TRUE)
      # Family (5): Methanosaetaceae → Methanotrichaceae — before Genus rename
      p[5L] <- gsub("Methanosaetaceae", "Methanotrichaceae", p[5L], fixed = TRUE)
      # Genus (6) + Species (7): Methanosaeta → Methanothrix (Oren 2014)
      p[6L] <- gsub("Methanosaeta",     "Methanothrix",      p[6L], fixed = TRUE)
      p[7L] <- gsub("Methanosaeta",     "Methanothrix",      p[7L], fixed = TRUE)
      # Order (4): Methanosarcinales → Methanotrichales, only when Family = Methanotrichaceae
      if (identical(p[5L], "Methanotrichaceae") && identical(p[4L], "Methanosarcinales"))
        p[4L] <- "Methanotrichales"
      paste(p, collapse = ";")
    }

    tmpl_lines <- readLines(YANG_TEMPLATE_FA)
    out_lines  <- character(0)
    cur_acc    <- NULL; cur_seq <- character(0)

    flush_yang <- function() {
      if (is.null(cur_acc)) return()
      taxid <- acc_to_taxid[[cur_acc]]
      lin   <- if (!is.null(taxid) && nzchar(taxid) && !is.na(taxid))
                 taxid_lineage[[taxid]] else NULL
      if (!is.null(lin)) {
        tax_str <- paste(vapply(names(LRANKS_MAP), function(rk) {
          val <- lin[[rk]]; if (!is.null(val) && nzchar(trimws(val))) trimws(val)
          else "Unclassified"
        }, character(1L)), collapse = ";")
      } else {
        orig <- yang_map[[cur_acc]]
        tax_str <- if (!is.null(orig) && nzchar(orig)) pad_to_7(orig)
                   else paste(rep("Unclassified", 7L), collapse = ";")
      }
      tax_str   <- normalise_yang_tax(tax_str)
      out_lines <<- c(out_lines, paste0(">", tax_str, ";", cur_acc),
                      paste(cur_seq, collapse = ""))
    }

    for (ln in tmpl_lines) {
      if (startsWith(ln, ">")) {
        flush_yang(); cur_acc <- trimws(substring(ln, 2L)); cur_seq <- character(0L)
      } else cur_seq <- c(cur_seq, ln)
    }
    flush_yang()

    writeLines(out_lines, YANG_TAXIZE_FA)
    db_log(sprintf("  Written: %s  (%d sequences)",
                   YANG_TAXIZE_FA, sum(startsWith(out_lines, ">"))))
  } else {
    db_log(sprintf("  %s already exists — skipping", YANG_TAXIZE_FA))
  }
} else {
  db_log("  STEP 9 skipped: mcrAtemplate.fasta or Yang_etal_2014_tax4mcrA.taxonomy.txt not found",
         "WARN")
}
