#!/usr/bin/env Rscript
# =============================================================================
# build_mcrA_db.R
# ===============
# Seed-and-Expand mcrA Reference Database Builder for DADA2  (R version)
#
# Mirrors the Python workflow in build_mcrA_db_final.py (v2).
#
# Workflow:
#   0.  Cluster seed template at 95 % identity  (cd-hit-est via system call)
#   1.  Load seed sequences
#   2.  BLAST seeds against NCBI nt  (rentrez / NCBI BLAST API)
#   3.  Fetch full sequences + GenBank records from NCBI  (rentrez)
#   4.  Fetch taxonomy; correct names with taxize against current NCBI records
#   5.  HMMER gene-coverage filter + trimming  (hmmsearch via system call)
#   6.  Merge BLAST pool with yang_2014_taxize.fasta and mcrA_ncbi_genome_db.fasta;
#         deduplicate exact-sequence duplicates by retaining the entry with the
#         most complete lineage; append accession to the Species field
#   7.  Alignment-based gap trimming  (MAFFT via system call)
#         – trim alignment columns with ≥ GAP_THRESHOLD gap frequency
#         – remove sequences ending / starting > TRIM_MIN_COVERAGE dense
#           columns from the block boundaries
#   8.  Write mcrA_ncbi_nt_db.fasta  (NCBI taxonomy)
#   9.  Phylogenetic taxonomy curation  →  mcrA_ncbi_nt_cur_db.fasta
#         9a. Label cleaning before the tree step (removes local disagreements
#             between near-identical references):
#             – collapse unknown placeholders ("<X>_unclassified", "uncultured…")
#               to "Unclassified", the only value step 9b treats as unknown
#             – vsearch all-vs-all neighbours; per-rank % identity cut-off =
#               lowest identity above which ≥ PRECISION_TARGET of original-record
#               neighbour pairs share the taxon (non-decreasing kingdom → genus)
#             – each name takes the parent carried by ≥ MIN_AGREEMENT of the
#               original records with that name
#             – a named label is replaced when ≥ MIN_SUPPORT_FIX neighbours above
#               the cut-off agree on another taxon with ≥ MIN_AGREEMENT support,
#               and blanked when no neighbour supports it; genome-derived records
#               (GCA_/GCF_) are never changed
#         9b. Tree-based curation (unchanged rules)
#             – FastTree GTR tree from accepted sequences
#             – Label each internal node at rank R with taxon X only if:
#                 (i)   all classified tips in the subtree agree unanimously on X
#                 (ii)  max intra-subtree patristic distance ≤ PHYLO_RANK_THRESHOLDS
#                 (iii) ALL tips in the full tree carrying X at rank R fall within
#                       the subtree (strict monophyly; paraphyletic nodes rejected)
#             – For each "Unclassified" tip: walk tip → root; assign ranks from
#               the deepest ancestor that carries a qualifying label; fill coarser
#               ranks still unresolved by unanimous consensus of known sequences
#               in that ancestor's subtree
#             – Taxize re-verification + higher-rank lineage fill
#         9c. After the tree step
#             – neighbour check (as 9a) of the labels added by 9b
#             – fill ranks still "Unclassified" when ≥ MIN_SUPPORT_FILL
#               original-record neighbours agree (≥ MIN_AGREEMENT) and the parent
#               matches the sequence's own parent label
#             – lineage consistency (as 9a)
#         9d. Write the curated FASTA (8 fields), the dada2 assignTaxonomy
#             training file (Kingdom…Genus), the dada2 addSpecies file and a
#             change log with the evidence for every label change
#   10. Write accession correspondence TSV
#
# Required R packages (install once):
#   install.packages(c("rentrez", "taxize", "Biostrings", "stringr",
#                      "dplyr", "readr", "httr"))
#   # httr is used for direct NCBI BLAST CGI API calls
#   # Biostrings is a Bioconductor package:
#   if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
#   BiocManager::install("Biostrings")
#
# External tools (must be on PATH or paths set below):
#   conda install -c bioconda hmmer mafft blast cd-hit fasttree vsearch
#
# Usage:
#   Rscript build_mcrA_db.R
# =============================================================================

suppressPackageStartupMessages({
  library(rentrez)          # NCBI Entrez API (efetch, esearch, etc.)
  library(taxize)           # taxonomy name resolution / correction
  library(Biostrings)       # DNA sequence manipulation
  library(stringr)          # string utilities
  library(dplyr)            # data-frame manipulation
  library(readr)            # write_tsv
  library(httr)             # HTTP for NCBI BLAST CGI API and efetch
  library(xml2)             # XML parsing for taxonomy records
  library(ape)              # phylogenetic tree I/O + patristic distances (Step 9)
  library(phangorn)         # Descendants / Ancestors traversal (Step 9)
})

# =============================================================================
# CONFIGURATION
# =============================================================================

# External databases merged before Step 6
# (built by build_mcrA_db_ncbi_annot.R — run that script first)
YANG_TAXIZE_FILE  <- "yang_2014_taxize.fasta"    # Yang 2014 with NCBI-enriched lineages
NCBI_ANNOT_FILE   <- "mcrA_ncbi_genome_db.fasta"  # NCBI genome annotation-based database

SEED_TEMPLATE    <- YANG_TAXIZE_FILE     # yang_2014_taxize.fasta replaces mcrAtemplate.fasta
SEED_FILE        <- "cdhit_out"
OUTPUT_FILE      <- "mcrA_ncbi_nt_db.fasta"
OUTPUT_CORRESPONDENCE <- "mcrA_accession_correspondence.tsv"
LOG_FILE         <- "build_mcrA_db_R.log"

ENTREZ_EMAIL     <- Sys.getenv("NCBI_EMAIL",    "your.email@example.com")  # required by NCBI
ENTREZ_API_KEY   <- Sys.getenv("NCBI_API_KEY", "")  # optional but raises rate limit to 10 req/s
if (nchar(ENTREZ_API_KEY) > 0 && !nzchar(Sys.getenv("ENTREZ_KEY")))
  Sys.setenv(ENTREZ_KEY = ENTREZ_API_KEY)   # taxize reads ENTREZ_KEY
# To set permanently, add to ~/.Renviron:
#   NCBI_EMAIL=your.email@example.com
#   NCBI_API_KEY=<your key from https://www.ncbi.nlm.nih.gov/account/>
# BLAST parameters (web BLAST against NCBI nt)
BLAST_DB          <- "nt"
MIN_PERC_IDENTITY <- 80
MAX_HITS_PER_SEED <- 500
ENTREZ_QUERY      <- "Archaea[Organism] AND mcrA[Gene Name]"

# GenBank CDS extraction keywords (fallback when HMM_PROFILE is NULL)
EXTRACT_COMPLETE_GENE <- TRUE
MCRA_GENE_KEYWORDS    <- c("mcra", "methyl-coenzyme m reductase",
                            "methyl coenzyme m reductase", "mcr alpha", "mcrα")

# Sequence quality filters
MIN_LENGTH        <- 300
MAX_LENGTH        <- 2000
MAX_AMBIGUOUS_PCT <- 0.05

# NCBI rate limiting (seconds between requests)
SLEEP_SEC         <- if (nchar(ENTREZ_API_KEY) > 0) 0.12 else 0.4
FETCH_BATCH_SIZE  <- 200   # smaller batches for R's rentrez

# ── HMMER ─────────────────────────────────────────────────────────────────────
HMM_PROFILE       <- "TIGR03256.hmm"   # set to NULL to skip
HMM_MIN_COVERAGE  <- 0.20
HMMSEARCH_PATH    <- "hmmsearch"
HMMER_THREADS     <- 24
HMMER_EVALUE      <- "0.00001"
HMMER_MAX_NT_LEN  <- 50000
HMMER_CACHE_FILE  <- "hmm_trimmed.rds" # Step 5 checkpoint — delete to re-run HMMER

# ── cd-hit-est (Step 0 seed clustering only) ──────────────────────────────────
CDHIT_PATH           <- "cd-hit-est"
CDHIT_SEED_IDENTITY  <- 0.95
CDHIT_SEED_WORD_SIZE <- 8
CDHIT_THREADS        <- 24

# ── BLAST result caching ──────────────────────────────────────────────────────
# Two-level cache to survive downstream failures without re-BLASTing:
#   Level 1 — per-seed XML:  blast_xml_cache/<seed_id>.xml
#             If a seed's XML is already on disk, its accessions are parsed
#             locally and no remote BLAST call is made for that seed.
#   Level 2 — merged accession list:  BLAST_CACHE_FILE
#             Written after ALL seeds finish.  On the next run the whole
#             BLAST step is skipped and accessions are read from this file.
# Delete BLAST_CACHE_FILE (or the whole blast_xml_cache/ dir) to force a
# fresh BLAST run.
BLAST_CACHE_FILE <- "blast_accessions_cache.txt"   # one accession per line
BLAST_CACHE_DIR  <- "blast_xml_cache"              # per-seed raw XML files

# ── GenBank sequence caching (Step 3) ────────────────────────────────────────
# Two-level cache analogous to the BLAST cache:
#   Level 1 — per-batch raw GenBank text:  genbank_batch_cache/<batch_idx>.gb
#             Written immediately after each successful efetch batch so that a
#             mid-run failure does not discard already-fetched data.
#   Level 2 — parsed seq_df as RDS:  genbank_seq_df.rds
#             Written after parse_genbank_records() completes.  On the next run
#             the whole fetch+parse step is skipped by loading this RDS.
# Delete genbank_seq_df.rds (or the whole genbank_batch_cache/ dir) to force a
# fresh fetch.
GENBANK_CACHE_DIR <- "genbank_batch_cache"  # per-batch raw GenBank text files
SEQ_CACHE_FILE    <- "genbank_seq_df.rds"   # parsed seq_df (survives downstream failures)
TAX_CACHE_FILE    <- "taxid_lineage.rds"    # resolved taxid → lineage list (Step 4)

# (Steps 5b/5c removed — taxonomy now comes from yang_2014_taxize.fasta and
#  mcrA_ncbi_genome_db.fasta merged in Step 6; TAX4MCRA_FILE, YANG_FASTA_FILE,
#  YANG_TAXONOMY_FILE, and their associated BLAST vars are no longer needed.)

# ── Taxonomy name correction via taxize ───────────────────────────────────────
CORRECT_TAXONOMY <- TRUE           # verify names against current NCBI taxonomy

# ── Alignment-based gap trimming (MAFFT) ─────────────────────────────────────
MAFFT_PATH        <- "mafft"
MAFFT_THREADS     <- 24
GAP_THRESHOLD     <- 0.50    # columns with ≥ 50 % gaps → sparse
TRIM_MIN_COVERAGE <- 50      # remove seqs missing > 50 dense columns at either end

# ── Phylogenetic taxonomy curation (Step 9) ────────────────────────────────────
# Builds a tree from the final sequences, labels internal nodes where all known
# descendants agree on a taxon AND the max patristic distance is below the
# rank threshold, then propagates labels to "Unclassified" tips within those
# monophyletic clades.  Writes both a curated FASTA and the tree (Newick) for
# external review (e.g., iTOL, FigTree).
#
# Thresholds are PATRISTIC distances — sums of GTR branch lengths along the path
# between two tips in the FastTree tree (ape::cophenetic.phylo).  They are NOT
# pairwise K80 sequence distances.  Patristic distances are typically 4–8× larger
# than equivalent K80 pairwise distances for the same evolutionary divergence, so
# values above 1.0 are normal at class / phylum level.
#
# The values below were calibrated empirically with calibrate_mcra_thresholds.R
# using the Youden-J method on mcrAtemplate.fasta + tax4mcrA.taxonomy.
# Species and Phylum are omitted: Species returns J = 0 (AUC ≈ 0.08, threshold
# −Inf) and Phylum returns J = 0.03 (AUC ≈ 0.50), indicating that neither rank
# boundary can be reliably inferred from patristic distance in this dataset.

PHYLO_CURATE_TAXONOMY <- TRUE          # set FALSE to skip Step 9
FASTTREE_PATH         <- "~/anaconda3/bin/FastTree"  # full path; auto-resolved at runtime
TREE_FILE             <- "mcrA_ncbi_nt_db.nwk"         # Newick tree — open in FigTree / iTOL
OUTPUT_FILE_CURATED   <- "mcrA_ncbi_nt_cur_db.fasta"    # taxonomy-curated DADA2 FASTA

# ── Step 9a / 9c: neighbour-guided label cleaning and filling ─────────────────
VSEARCH_PATH       <- "vsearch"
VSEARCH_THREADS    <- 24
NEIGHBOUR_MIN_ID   <- 0.80   # lowest identity reported by vsearch
NEIGHBOUR_MAX_HITS <- 100    # --maxaccepts
MIN_AGREEMENT      <- 0.90   # fraction of voting neighbours that must agree
MIN_SUPPORT_FIX    <- 3      # neighbours needed to replace a named label
MIN_SUPPORT_FILL   <- 2      # neighbours needed to fill an "Unclassified" rank
CALIB_MIN_PAIRS    <- 30     # positive and negative pairs needed per rank
PRECISION_TARGET   <- 0.98   # share of neighbour pairs above the cut-off that
                             # must carry the same taxon (in mcrA_ncbi_nt_db
                             # 1-2 % of original-record pairs disagree even at
                             # 100 % identity, so 0.99 is not reachable at
                             # order-genus level)
FAST_NODE_LABELLING <- TRUE  # exact re-implementations of node labelling and
                             # leaf-to-root propagation (the original functions
                             # need a 27,942 x 27,942 patristic matrix and tip
                             # lists for every node); FALSE uses the originals
TRAIN_MIN_RANK     <- "phylum"  # dada2 training set: references unknown at this
                                # rank are left out
REFINE_UNRES_PATTERN <- "^(Unclassified|unclassified|uncultured|NA|)$|^uncultured|_unclassified$"
NEIGHBOURS_FILE              <- "mcrA_ncbi_nt_db_neighbours.tsv"
OUTPUT_FILE_CURATED_TRAIN    <- "mcrA_ncbi_nt_cur_db_dada2_toGenus.fasta"  # assignTaxonomy
OUTPUT_FILE_CURATED_SPECIES  <- "mcrA_ncbi_nt_cur_db_dada2_species.fasta"  # addSpecies
OUTPUT_FILE_CURATION_CHANGES <- "mcrA_ncbi_nt_cur_db_changes.tsv"
OUTPUT_FILE_CURATION_CUTOFFS <- "mcrA_ncbi_nt_cur_db_identity_cutoffs.tsv"
OUTPUT_FILE_CURATION_SUMMARY <- "mcrA_ncbi_nt_cur_db_curation_summary.tsv"

PHYLO_RANK_THRESHOLDS <- c(           # max patristic distance per rank (GTR subst/site via FastTree)
  # species omitted — J=0.000, AUC=0.078: cannot discriminate by distance
  genus   = 0.3595,  # J=0.587  AUC=0.817  95%CI [0.3483, 0.3603]
  family  = 0.8045,  # J=0.280  AUC=0.662  95%CI [0.8028, 0.8077]  (use with caution)
  order   = 0.8135,  # J=0.610  AUC=0.803  95%CI [0.8005, 0.8139]
  class   = 1.1446   # J=0.439  AUC=0.778  95%CI [1.1356, 1.1481]
  # phylum omitted — J=0.032, AUC=0.504: effectively random discrimination
)

# DADA2 taxonomy levels (accession appended separately after Species)
TAX_LEVELS <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")

# NCBI rank → DADA2 level
NCBI_RANK_MAP <- c(
  superkingdom = "kingdom",
  domain       = "kingdom",
  phylum       = "phylum",
  class        = "class",
  order        = "order",
  family       = "family",
  genus        = "genus",
  species      = "species"
)

KNOWN_KINGDOMS <- c(
  archaea   = "Archaea",
  bacteria  = "Bacteria",
  eukaryota = "Eukaryota",
  viruses   = "Viruses"
)

# =============================================================================
# LOGGING HELPER
# =============================================================================

log_con <- file(LOG_FILE, open = "w")

mcra_log <- function(msg, level = "INFO") {
  ts  <- format(Sys.time(), "%H:%M:%S")
  txt <- sprintf("%s  %-8s  %s", ts, level, msg)
  message(txt)
  cat(txt, "\n", file = log_con, append = TRUE)
}

on.exit({
  try(close(log_con), silent = TRUE)
}, add = TRUE)

# =============================================================================
# FASTA HELPERS
# =============================================================================

read_fasta <- function(path) {
  # Returns a named character vector: name → sequence
  lines <- readLines(path)
  header_idx <- which(startsWith(lines, ">"))
  seqs  <- character(length(header_idx))
  names <- sub("^>", "", lines[header_idx])
  for (i in seq_along(header_idx)) {
    start <- header_idx[i] + 1L
    end   <- if (i < length(header_idx)) header_idx[i + 1L] - 1L else length(lines)
    seqs[i] <- paste(lines[start:end], collapse = "")
  }
  names(seqs) <- names
  seqs
}

write_fasta <- function(seqs, path, headers = NULL) {
  # seqs: named character vector; headers overrides names if provided
  ids <- if (!is.null(headers)) headers else names(seqs)
  lines <- character(length(seqs) * 2)
  for (i in seq_along(seqs)) {
    hdr <- ids[i]
    if (!startsWith(hdr, ">")) hdr <- paste0(">", hdr)
    lines[2 * i - 1] <- hdr
    lines[2 * i]     <- seqs[i]
  }
  writeLines(lines, path)
}


# =============================================================================
# STEP 0 — Cluster seed template at 95 % with cd-hit-est
# =============================================================================

cluster_seeds <- function(template_file, output_file) {
  if (!file.exists(template_file)) {
    stop(sprintf("Seed template not found: %s", template_file))
  }
  cmd <- sprintf(
    "%s -i %s -o %s -c %.2f -n %d -T %d -M 0 -d 0",
    CDHIT_PATH, shQuote(template_file), shQuote(output_file),
    CDHIT_SEED_IDENTITY, CDHIT_SEED_WORD_SIZE, CDHIT_THREADS
  )
  mcra_log(sprintf("  Command: %s", cmd))
  ret <- system(cmd, ignore.stdout = TRUE, ignore.stderr = FALSE)
  if (ret != 0) stop("cd-hit-est failed")
  n_rep <- sum(startsWith(readLines(output_file), ">"))
  mcra_log(sprintf("  Seed clustering: %d representatives → %s", n_rep, output_file))
  invisible(TRUE)
}


# =============================================================================
# STEP 1 — Load seeds
# =============================================================================

load_seeds <- function(fasta_file) {
  if (!file.exists(fasta_file)) stop(sprintf("Seed file not found: %s", fasta_file))
  seeds <- read_fasta(fasta_file)
  if (length(seeds) == 0) stop(sprintf("No sequences in %s", fasta_file))
  mcra_log(sprintf("Loaded %d seed sequences from '%s'", length(seeds), fasta_file))
  seeds
}


# =============================================================================
# STEP 2 — BLAST seeds against NCBI nt  (BLAST CGI API via httr)
# =============================================================================
#
# rentrez::entrez_post() is for the Entrez History Server, NOT for BLAST.
# NCBI web BLAST requires direct HTTP calls to the BLAST CGI endpoint.
#
BLAST_CGI <- "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"

blast_seed_ncbi <- function(seed_seq, seed_id, xml_cache_path = NULL) {
  mcra_log(sprintf("  BLASTing '%s' (%d bp) ...", seed_id, nchar(seed_seq)))

  # ── 1. Submit job ──────────────────────────────────────────────────────────
  put_params <- list(
    CMD          = "Put",
    PROGRAM      = "blastn",
    DATABASE     = BLAST_DB,
    QUERY        = as.character(seed_seq),
    HITLIST_SIZE = as.character(MAX_HITS_PER_SEED),
    PERC_IDENT   = as.character(as.integer(MIN_PERC_IDENTITY)),
    FORMAT_TYPE  = "XML"
  )
  if (nchar(ENTREZ_QUERY) > 0) put_params$ENTREZ_QUERY <- ENTREZ_QUERY

  put_resp <- tryCatch(
    httr::POST(BLAST_CGI, body = put_params, encode = "form",
               httr::timeout(120)),
    error = function(e) {
      mcra_log(sprintf("  BLAST submit failed: %s", e$message), "WARN")
      NULL
    }
  )
  if (is.null(put_resp) || httr::http_error(put_resp)) {
    mcra_log("  BLAST submit: HTTP error", "WARN")
    return(character(0))
  }

  # Extract RID from the response HTML
  # NCBI embeds: <!--QBlastInfoBegin  RID = XXXXXXXX  ...  QBlastInfoEnd-->
  put_txt <- httr::content(put_resp, as = "text", encoding = "UTF-8")
  rid_m   <- regmatches(put_txt,
                         regexpr("RID\\s*=\\s*(\\S+)", put_txt, perl = TRUE))
  if (length(rid_m) == 0 || !nzchar(rid_m)) {
    mcra_log("  BLAST submit: could not parse RID from response", "WARN")
    return(character(0))
  }
  rid <- trimws(sub("RID\\s*=\\s*", "", rid_m))
  mcra_log(sprintf("  BLAST RID: %s — polling for results...", rid))

  # ── 2. Poll until ready (max ~13 min) ────────────────────────────────────
  Sys.sleep(15)
  ready <- FALSE
  for (attempt in seq_len(80)) {
    Sys.sleep(10)
    status_resp <- tryCatch(
      httr::GET(BLAST_CGI,
                query = list(CMD = "Get", RID = rid,
                             FORMAT_OBJECT = "SearchInfo"),
                httr::timeout(45),
                httr::config(ssl_verifypeer = 0L, ssl_verifyhost = 0L)),
      error = function(e) NULL
    )
    if (!is.null(status_resp)) {
      status_txt <- httr::content(status_resp, as = "text", encoding = "UTF-8")
      if (grepl("Status=READY",   status_txt, fixed = TRUE)) { ready <- TRUE; break }
      if (grepl("Status=FAILED",  status_txt, fixed = TRUE) ||
          grepl("Status=UNKNOWN", status_txt, fixed = TRUE)) {
        mcra_log(sprintf("  BLAST job %s failed/expired", rid), "WARN")
        return(character(0))
      }
    }
    if (attempt %% 6 == 0)
      mcra_log(sprintf("  Still waiting for BLAST RID %s (attempt %d)...", rid, attempt))
  }
  if (!ready) {
    mcra_log(sprintf("  BLAST job %s timed out after polling", rid), "WARN")
    return(character(0))
  }

  # ── 3. Retrieve XML results (stream to file; retry on TLS/network errors) ─
  # Large BLAST result sets can exceed 10 MB.  Streaming to disk with
  # httr::write_disk() avoids in-memory buffering and is more robust against
  # GnuTLS non-proper-termination errors (-110) on Linux.
  tmp_xml    <- tempfile(fileext = ".blast.xml")
  xml_str    <- NULL
  max_tries  <- 4

  for (try_n in seq_len(max_tries)) {
    if (try_n > 1) {
      wait_sec <- 15 * try_n
      mcra_log(sprintf("  BLAST retrieve retry %d/%d for RID %s (waiting %ds)...",
                       try_n, max_tries, rid, wait_sec))
      Sys.sleep(wait_sec)
    }

    result_resp <- tryCatch(
      httr::GET(
        BLAST_CGI,
        query  = list(CMD = "Get", RID = rid, FORMAT_TYPE = "XML"),
        httr::write_disk(tmp_xml, overwrite = TRUE),
        httr::timeout(300),                    # 5-min transfer window
        httr::config(ssl_verifypeer = 0L,      # tolerate non-clean TLS shutdown
                     ssl_verifyhost = 0L)
      ),
      error = function(e) {
        mcra_log(sprintf("  BLAST retrieve error (try %d): %s", try_n, e$message), "WARN")
        NULL
      }
    )

    if (!is.null(result_resp) && !httr::http_error(result_resp) &&
        file.exists(tmp_xml) && file.info(tmp_xml)$size > 0) {
      xml_str <- paste(readLines(tmp_xml, warn = FALSE), collapse = "\n")
      break
    }

    # If the file was written but contains an error page, report it
    if (file.exists(tmp_xml) && file.info(tmp_xml)$size > 0) {
      head_txt <- paste(readLines(tmp_xml, n = 5, warn = FALSE), collapse = " ")
      mcra_log(sprintf("  BLAST retrieve (try %d) bad response: %s", try_n,
                       substr(head_txt, 1, 200)), "WARN")
    }
  }

  if (is.null(xml_str) || !nzchar(xml_str)) {
    mcra_log(sprintf("  BLAST retrieve failed for RID %s after %d tries", rid, max_tries), "WARN")
    return(character(0))
  }

  # ── 4. Persist XML to per-seed cache ─────────────────────────────────────
  if (!is.null(xml_cache_path)) {
    tryCatch({
      writeLines(readLines(tmp_xml, warn = FALSE), xml_cache_path)
      mcra_log(sprintf("  Cached XML → %s", xml_cache_path))
    }, error = function(e) {
      mcra_log(sprintf("  Could not write XML cache: %s", e$message), "WARN")
    })
  }

  # ── 5. Parse accessions from XML ─────────────────────────────────────────
  acc_hits    <- str_extract_all(xml_str,
                                  "(?<=<Hit_accession>)[^<]+(?=</Hit_accession>)")[[1]]
  unique_accs <- unique(str_trim(acc_hits))
  mcra_log(sprintf("  → %d unique accessions", length(unique_accs)))
  unique_accs
}

blast_all_seeds <- function(seeds) {

  # ── Level-2 cache: full merged accession list ─────────────────────────────
  if (file.exists(BLAST_CACHE_FILE)) {
    cached <- readLines(BLAST_CACHE_FILE, warn = FALSE)
    cached <- cached[nzchar(trimws(cached))]
    mcra_log(sprintf(
      "  Found merged BLAST cache '%s' — loading %d accessions (skipping remote BLAST).",
      BLAST_CACHE_FILE, length(cached)))
    mcra_log("  Delete this file to force a fresh BLAST run.")
    return(cached)
  }

  # ── Level-1 cache: per-seed XML directory ────────────────────────────────
  if (!dir.exists(BLAST_CACHE_DIR)) {
    dir.create(BLAST_CACHE_DIR, recursive = TRUE)
    mcra_log(sprintf("  Created per-seed XML cache dir: %s/", BLAST_CACHE_DIR))
  }

  # Helper: safe filename for a seed ID
  safe_seed_id <- function(sid)
    gsub("[^A-Za-z0-9._-]", "_", sid, perl = TRUE)

  parse_accs_from_xml <- function(xml_path) {
    xml_str  <- paste(readLines(xml_path, warn = FALSE), collapse = "\n")
    acc_hits <- str_extract_all(xml_str,
                                 "(?<=<Hit_accession>)[^<]+(?=</Hit_accession>)")[[1]]
    unique(str_trim(acc_hits))
  }

  all_accs <- character(0)

  for (i in seq_along(seeds)) {
    seed_id      <- names(seeds)[i]
    xml_cache    <- file.path(BLAST_CACHE_DIR,
                              paste0(safe_seed_id(seed_id), ".xml"))
    mcra_log(sprintf("[BLAST %d/%d]", i, length(seeds)))

    # Check per-seed XML cache before going remote
    if (file.exists(xml_cache) && file.info(xml_cache)$size > 100) {
      mcra_log(sprintf("  Using cached XML for '%s'", seed_id))
      accs <- tryCatch(parse_accs_from_xml(xml_cache),
                       error = function(e) {
                         mcra_log(sprintf("  Cache parse failed (%s) — re-BLASTing", e$message), "WARN")
                         blast_seed_ncbi(seeds[i], seed_id, xml_cache)
                       })
    } else {
      accs <- blast_seed_ncbi(seeds[i], seed_id, xml_cache)
    }

    all_accs <- unique(c(all_accs, accs))
    mcra_log(sprintf("  Cumulative unique accessions: %d", length(all_accs)))
    Sys.sleep(SLEEP_SEC)
  }

  # ── Save merged accession cache (Level 2) ────────────────────────────────
  tryCatch({
    writeLines(all_accs, BLAST_CACHE_FILE)
    mcra_log(sprintf("  Saved %d accessions to merged cache '%s'",
                     length(all_accs), BLAST_CACHE_FILE))
  }, error = function(e) {
    mcra_log(sprintf("  Could not write merged cache: %s", e$message), "WARN")
  })

  mcra_log(sprintf("BLAST complete. Total unique accessions: %d", length(all_accs)))
  all_accs
}


# =============================================================================
# STEP 3 — Fetch GenBank records from NCBI  (rentrez)
# =============================================================================

fetch_sequences <- function(accessions) {
  # Fetch GenBank flat-file records from NCBI in batches.
  # Uses httr::GET directly against efetch for TLS robustness (ssl_verifypeer=0),
  # and writes each batch to a per-batch cache file in GENBANK_CACHE_DIR so that
  # a mid-run failure does not throw away already-downloaded data.
  #
  # Returns a list of lists, each with element $gb_text (character string).

  total  <- length(accessions)
  records <- list()
  mcra_log(sprintf("Fetching %d sequences from NCBI (batch=%d)...", total, FETCH_BATCH_SIZE))

  # Ensure per-batch cache directory exists
  if (!dir.exists(GENBANK_CACHE_DIR)) {
    dir.create(GENBANK_CACHE_DIR, recursive = TRUE)
    mcra_log(sprintf("  Created GenBank batch cache dir: %s/", GENBANK_CACHE_DIR))
  }

  efetch_base <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

  batch_idx <- 0L
  for (i in seq(1, total, by = FETCH_BATCH_SIZE)) {
    batch_idx <- batch_idx + 1L
    end_i     <- min(i + FETCH_BATCH_SIZE - 1L, total)
    batch     <- accessions[i:end_i]

    cache_path <- file.path(GENBANK_CACHE_DIR,
                            sprintf("batch_%05d.gb", batch_idx))

    # ── Use per-batch cache if available ──────────────────────────────────────
    if (file.exists(cache_path) && file.info(cache_path)$size > 0L) {
      mcra_log(sprintf("  Batch %d/%d (%d seqs) — loading from cache",
                       batch_idx,
                       ceiling(total / FETCH_BATCH_SIZE), length(batch)))
      gb_txt <- paste(readLines(cache_path, warn = FALSE), collapse = "\n")
      records[[length(records) + 1L]] <- list(gb_text = gb_txt)
      next
    }

    mcra_log(sprintf("  Fetching sequences %d–%d of %d (batch %d)...",
                     i, end_i, total, batch_idx))

    # ── httr GET to efetch (retry up to 3 times) ──────────────────────────────
    query_params <- list(
      db      = "nuccore",
      id      = paste(batch, collapse = ","),
      rettype = "gb",
      retmode = "text"
    )
    if (nchar(ENTREZ_API_KEY) > 0L)
      query_params$api_key <- ENTREZ_API_KEY

    gb_txt <- NULL
    for (try_n in seq_len(3L)) {
      if (try_n > 1L) {
        wait_sec <- 15L * try_n
        mcra_log(sprintf("  Retry %d/3 for batch %d (waiting %ds)...",
                         try_n, batch_idx, wait_sec), "WARN")
        Sys.sleep(wait_sec)
      }

      resp <- tryCatch(
        httr::GET(
          efetch_base,
          query  = query_params,
          httr::timeout(300L),
          httr::config(ssl_verifypeer = 0L, ssl_verifyhost = 0L)
        ),
        error = function(e) {
          mcra_log(sprintf("  Batch %d fetch error (try %d): %s",
                           batch_idx, try_n, e$message), "WARN")
          NULL
        }
      )

      if (!is.null(resp) && !httr::http_error(resp)) {
        gb_txt <- httr::content(resp, as = "text", encoding = "UTF-8")
        if (!is.null(gb_txt) && nchar(gb_txt) > 0L) break
      } else if (!is.null(resp)) {
        mcra_log(sprintf("  Batch %d HTTP error %d (try %d)",
                         batch_idx, httr::status_code(resp), try_n), "WARN")
      }
    }

    if (is.null(gb_txt) || nchar(gb_txt) == 0L) {
      mcra_log(sprintf("  Batch %d fetch failed after retries — skipping", batch_idx), "WARN")
      Sys.sleep(SLEEP_SEC)
      next
    }

    # ── Persist to per-batch cache ────────────────────────────────────────────
    tryCatch(
      writeLines(gb_txt, cache_path),
      error = function(e)
        mcra_log(sprintf("  Could not write batch cache %s: %s", cache_path, e$message), "WARN")
    )

    records[[length(records) + 1L]] <- list(gb_text = gb_txt)
    Sys.sleep(SLEEP_SEC)
  }

  mcra_log(sprintf("Fetch complete. %d/%d batches retrieved.",
                   length(records), batch_idx))
  records   # list of lists; each has $gb_text (character string)
}


# =============================================================================
# GenBank text parser helpers
# =============================================================================

parse_genbank_records <- function(gb_text_list) {
  # Returns a data.frame with columns: accession, taxid, sequence, description
  rows <- list()

  for (gb_text in gb_text_list) {
    # Split into individual records on "//" boundary
    record_blocks <- strsplit(gb_text$gb_text, "\n//\n?")[[1]]
    for (blk in record_blocks) {
      if (!nchar(trimws(blk))) next

      # Accession
      acc_m <- regmatches(blk, regexpr("ACCESSION\\s+(\\S+)", blk, perl = TRUE))
      if (length(acc_m) == 0) next
      accession <- trimws(sub("ACCESSION\\s+", "", acc_m))

      # Taxid from /db_xref="taxon:XXXXXX"
      taxid_m <- regmatches(blk, regexpr('/db_xref="taxon:(\\d+)"', blk, perl = TRUE))
      taxid   <- if (length(taxid_m) > 0)
                   sub('/db_xref="taxon:(\\d+)"', "\\1", taxid_m, perl = TRUE)
                 else NA_character_

      # Sequence (ORIGIN section)
      # NOTE: regexpr("ORIGIN.*") only matches to end-of-line because "." in R
      # does not cross newlines by default.  Use string-position arithmetic instead.
      seq_str    <- ""
      origin_pos <- regexpr("ORIGIN", blk, fixed = TRUE)
      if (origin_pos > 0L) {
        after_origin <- substring(blk, origin_pos[1L] + 6L)   # skip "ORIGIN"
        nl_pos       <- regexpr("\n", after_origin, fixed = TRUE)
        seq_raw      <- if (nl_pos > 0L) substring(after_origin, nl_pos + 1L) else after_origin
        seq_str      <- toupper(gsub("[^A-Za-z]", "", seq_raw))
      }
      if (nchar(seq_str) == 0L) next

      # mcrA CDS extraction (if requested)
      cds_seq <- NA_character_
      if (EXTRACT_COMPLETE_GENE) {
        cds_sec <- regmatches(blk, regexpr("(?s)\\s+CDS\\s+[^\n]+(?:\\n\\s+[^A-Z][^\n]+)*", blk, perl = TRUE))
        if (length(cds_sec) > 0) {
          blk_low <- tolower(cds_sec)
          if (any(sapply(MCRA_GENE_KEYWORDS, function(kw) grepl(kw, blk_low, fixed = TRUE)))) {
            # Try to extract coordinates  (simplified: use full seq_str for now)
            # A production implementation would parse complement(join(..)) coordinates.
            cds_seq <- seq_str
          }
        }
      }

      rows[[length(rows) + 1]] <- data.frame(
        accession   = accession,
        taxid       = taxid,
        sequence    = seq_str,
        cds_seq     = if (!is.na(cds_seq)) cds_seq else seq_str,
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(rows) == 0) return(data.frame())
  bind_rows(rows)
}


# =============================================================================
# STEP 4 — Fetch and correct taxonomy
# =============================================================================

parse_ncbi_taxonomy <- function(tax_record_xml) {
  # Extract lineage from rentrez XML list-of-lists structure
  lineage <- setNames(rep("Unclassified", length(TAX_LEVELS)), TAX_LEVELS)

  if (is.null(tax_record_xml)) return(lineage)

  # LineageEx contains ranked taxa
  lin_ex <- tax_record_xml$LineageEx
  if (!is.null(lin_ex)) {
    for (node in lin_ex) {
      rank <- tolower(node$Rank)
      if (rank %in% names(NCBI_RANK_MAP)) {
        dada2_rank <- NCBI_RANK_MAP[rank]
        lineage[dada2_rank] <- node$ScientificName
      }
    }
  }

  # Fill species if record itself is species-rank
  if (!is.null(tax_record_xml$Rank) && tolower(tax_record_xml$Rank) == "species") {
    lineage["species"] <- tax_record_xml$ScientificName
  }

  # Kingdom fallback from flat Lineage string
  if (lineage["kingdom"] == "Unclassified" && !is.null(tax_record_xml$Lineage)) {
    tokens <- tolower(str_trim(strsplit(tax_record_xml$Lineage, ";")[[1]]))
    for (tok in tokens) {
      if (tok %in% names(KNOWN_KINGDOMS)) {
        lineage["kingdom"] <- KNOWN_KINGDOMS[tok]
        break
      }
    }
  }

  lineage
}

fetch_taxonomy_batch <- function(taxids) {
  # Fetch NCBI taxonomy XML via httr::GET (TLS-hardened) and parse with xml2.
  # rentrez::entrez_fetch(..., parsed=TRUE) returns an XML externalptr that
  # cannot be subset with $ — this implementation avoids that entirely.

  unique_taxids <- unique(taxids[!is.na(taxids) & nzchar(taxids)])
  total         <- length(unique_taxids)
  taxid_lineage <- list()

  mcra_log(sprintf("Fetching taxonomy for %d unique taxon IDs...", total))
  if (total == 0L) return(taxid_lineage)

  efetch_base <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

  for (i in seq(1L, total, by = FETCH_BATCH_SIZE)) {
    end_i <- min(i + FETCH_BATCH_SIZE - 1L, total)
    batch <- unique_taxids[i:end_i]
    mcra_log(sprintf("  Fetching taxonomy %d–%d of %d...", i, end_i, total))

    tryCatch({
      query_params <- list(
        db      = "taxonomy",
        id      = paste(batch, collapse = ","),
        rettype = "xml",
        retmode = "xml"
      )
      if (nchar(ENTREZ_API_KEY) > 0L)
        query_params$api_key <- ENTREZ_API_KEY

      resp <- httr::GET(
        efetch_base,
        query  = query_params,
        httr::timeout(120L),
        httr::config(ssl_verifypeer = 0L, ssl_verifyhost = 0L)
      )

      if (httr::http_error(resp)) {
        mcra_log(sprintf("  Taxonomy batch HTTP %d — skipping", httr::status_code(resp)), "WARN")
        Sys.sleep(SLEEP_SEC)
        next
      }

      xml_txt <- httr::content(resp, as = "text", encoding = "UTF-8")
      doc     <- xml2::read_xml(xml_txt)

      # Top-level <Taxon> nodes directly under <TaxaSet>
      taxon_nodes <- xml2::xml_find_all(doc, "/TaxaSet/Taxon")

      for (node in taxon_nodes) {
        tid_node <- xml2::xml_find_first(node, "TaxId")
        if (inherits(tid_node, "xml_missing")) next
        tid <- xml2::xml_text(tid_node)
        if (is.na(tid) || !nzchar(tid)) next

        lineage <- setNames(rep("Unclassified", length(TAX_LEVELS)), TAX_LEVELS)

        # ── LineageEx: ranked ancestors ──────────────────────────────────────
        lin_nodes <- xml2::xml_find_all(node, "LineageEx/Taxon")
        for (ln in lin_nodes) {
          rank_n <- xml2::xml_find_first(ln, "Rank")
          name_n <- xml2::xml_find_first(ln, "ScientificName")
          if (inherits(rank_n, "xml_missing") || inherits(name_n, "xml_missing")) next
          rank_lc <- tolower(xml2::xml_text(rank_n))
          sci_nm  <- xml2::xml_text(name_n)
          if (rank_lc %in% names(NCBI_RANK_MAP) && nzchar(sci_nm))
            lineage[NCBI_RANK_MAP[rank_lc]] <- sci_nm
        }

        # ── Species from the record itself if species-rank ───────────────────
        rank_self_n <- xml2::xml_find_first(node, "Rank")
        name_self_n <- xml2::xml_find_first(node, "ScientificName")
        if (!inherits(rank_self_n, "xml_missing") &&
            !inherits(name_self_n, "xml_missing")) {
          if (tolower(xml2::xml_text(rank_self_n)) == "species")
            lineage["species"] <- xml2::xml_text(name_self_n)
        }

        # ── Kingdom fallback from flat <Lineage> string ──────────────────────
        if (lineage["kingdom"] == "Unclassified") {
          lin_str_n <- xml2::xml_find_first(node, "Lineage")
          if (!inherits(lin_str_n, "xml_missing")) {
            tokens <- tolower(str_trim(strsplit(xml2::xml_text(lin_str_n), ";")[[1]]))
            for (tok in tokens) {
              if (tok %in% names(KNOWN_KINGDOMS)) {
                lineage["kingdom"] <- KNOWN_KINGDOMS[tok]
                break
              }
            }
          }
        }

        taxid_lineage[[tid]] <- lineage
      }

      mcra_log(sprintf("  Parsed %d taxonomy records", length(taxon_nodes)))
    }, error = function(e) {
      mcra_log(sprintf("  Taxonomy batch fetch failed: %s", e$message), "WARN")
    })

    Sys.sleep(SLEEP_SEC)
  }

  mcra_log(sprintf("Taxonomy retrieval complete. %d taxon IDs resolved.", length(taxid_lineage)))
  taxid_lineage
}


# ── Taxonomy name correction with taxize ─────────────────────────────────────

.tax_correction_cache <- new.env(parent = emptyenv())

should_skip_correction <- function(name) {
  if (is.na(name) || name == "Unclassified") return(TRUE)
  name_low <- tolower(name)
  any(startsWith(name_low, c("uncultured", "unclassified", "unknown",
                              "environmental", "candidate", "metagenome")))
}

correct_lineage_names <- function(lineage) {
  if (!CORRECT_TAXONOMY) return(lineage)

  species <- lineage["species"]
  if (should_skip_correction(species)) return(lineage)

  # Check cache
  cache_key <- as.character(species)
  if (exists(cache_key, envir = .tax_correction_cache)) {
    cached <- get(cache_key, envir = .tax_correction_cache)
    return(if (!is.null(cached)) cached else lineage)
  }

  corrected <- tryCatch({
    # Use taxize::get_uid to look up the species in NCBI taxonomy
    uid_result <- taxize::get_uid(species, db = "ncbi",
                                  ask = FALSE, messages = FALSE)
    uid <- uid_result[[1]]

    if (is.na(uid)) {
      # Try name resolution across databases to find synonyms
      resolved <- tryCatch(
        taxize::gnr_resolve(names = species, data_source_ids = 4,  # NCBI
                            canonical = TRUE, best_match_only = TRUE,
                            fields = "all"),
        error = function(e) NULL
      )
      if (!is.null(resolved) && nrow(resolved) > 0) {
        matched_name <- resolved$matched_name2[1]
        uid_result2  <- taxize::get_uid(matched_name, db = "ncbi",
                                        ask = FALSE, messages = FALSE)
        uid <- uid_result2[[1]]
      }
    }

    if (is.na(uid)) {
      assign(cache_key, NULL, envir = .tax_correction_cache)
      return(lineage)
    }

    # Fetch full classification
    cls <- taxize::classification(uid, db = "ncbi")[[1]]
    if (is.null(cls) || nrow(cls) == 0) {
      assign(cache_key, NULL, envir = .tax_correction_cache)
      return(lineage)
    }

    # Build corrected lineage
    new_lineage <- setNames(rep("Unclassified", length(TAX_LEVELS)), TAX_LEVELS)
    rank_map_lower <- tolower(names(NCBI_RANK_MAP))

    for (j in seq_len(nrow(cls))) {
      rank_low <- tolower(cls$rank[j])
      if (rank_low %in% rank_map_lower) {
        dada2_rank <- NCBI_RANK_MAP[match(rank_low, rank_map_lower)]
        new_lineage[dada2_rank] <- cls$name[j]
      }
    }

    current_name <- cls$name[cls$rank == "species"]
    if (length(current_name) > 0 && current_name != species) {
      mcra_log(sprintf("  Taxonomy correction: '%s' → '%s'", species, current_name))
    }

    assign(cache_key, new_lineage, envir = .tax_correction_cache)
    new_lineage

  }, error = function(e) {
    assign(cache_key, NULL, envir = .tax_correction_cache)
    lineage
  })

  corrected
}


# ── Fill Unclassified higher ranks from the most specific known rank ─────────
#
# If a sequence is annotated at, say, Genus level but has Unclassified at
# Class/Order/Family, this function queries NCBI Taxonomy via taxize using the
# most specific known rank name and fills in the missing higher ranks.
#
# Consistency check: if any EXISTING non-Unclassified rank in the lineage
# conflicts with what taxize returns (e.g. Order = "Methanomicrobiales" but
# taxize reports "Methanobacteriales"), the fill is skipped and a warning is
# logged — existing annotations are never silently overridden.
#
fill_lineage_from_taxize <- function(lineage) {
  if (!CORRECT_TAXONOMY) return(lineage)

  # Find the most specific (finest) non-Unclassified rank
  known_idx <- which(TAX_LEVELS %in%
    TAX_LEVELS[!is.na(lineage[TAX_LEVELS]) & lineage[TAX_LEVELS] != "Unclassified"])
  if (length(known_idx) == 0L) return(lineage)

  finest_idx   <- max(known_idx)                  # last = finest rank
  query_rank   <- TAX_LEVELS[finest_idx]
  query_name   <- lineage[query_rank]

  if (should_skip_correction(query_name)) return(lineage)

  # Identify Unclassified ranks that are coarser than the finest known rank
  coarser_unclass <- TAX_LEVELS[seq_len(finest_idx - 1L)]
  coarser_unclass <- coarser_unclass[
    !is.na(lineage[coarser_unclass]) & lineage[coarser_unclass] == "Unclassified"
  ]
  if (length(coarser_unclass) == 0L) return(lineage)   # nothing to fill

  # Cache key scoped to rank + name to avoid collisions with correct_lineage_names
  cache_key <- paste0("fill__", query_rank, "__", query_name)
  if (exists(cache_key, envir = .tax_correction_cache)) {
    cached <- get(cache_key, envir = .tax_correction_cache)
    if (!is.null(cached)) {
      for (r in coarser_unclass)
        if (!is.na(cached[r]) && cached[r] != "Unclassified") lineage[r] <- cached[r]
    }
    return(lineage)
  }

  # Query NCBI via taxize
  taxize_lin <- tryCatch({
    uid_res <- taxize::get_uid(query_name, db = "ncbi",
                               ask = FALSE, messages = FALSE)
    uid <- uid_res[[1L]]
    if (is.na(uid)) stop("no uid")
    cls <- taxize::classification(uid, db = "ncbi")[[1L]]
    if (is.null(cls) || nrow(cls) == 0L) stop("empty classification")

    new_lin      <- setNames(rep("Unclassified", length(TAX_LEVELS)), TAX_LEVELS)
    rl           <- tolower(names(NCBI_RANK_MAP))
    for (j in seq_len(nrow(cls))) {
      rk_low <- tolower(cls$rank[j])
      if (rk_low %in% rl)
        new_lin[NCBI_RANK_MAP[match(rk_low, rl)]] <- cls$name[j]
    }
    new_lin
  }, error = function(e) {
    assign(cache_key, NULL, envir = .tax_correction_cache)
    NULL
  })

  if (is.null(taxize_lin)) return(lineage)
  assign(cache_key, taxize_lin, envir = .tax_correction_cache)

  # Consistency check: taxize must agree with every existing non-Unclassified rank
  for (r in TAX_LEVELS) {
    exist <- lineage[r]; tval <- taxize_lin[r]
    if (!is.na(exist) && exist != "Unclassified" &&
        !is.na(tval)  && tval  != "Unclassified" &&
        tolower(trimws(exist)) != tolower(trimws(tval))) {
      mcra_log(
        sprintf("  [fill_lineage] %s=%s: taxize conflict at %s ('%s' vs '%s') — fill skipped",
                query_rank, query_name, r, exist, tval), "WARN")
      return(lineage)   # abort entire fill to avoid partial inconsistency
    }
  }

  # Fill Unclassified coarser ranks
  for (r in coarser_unclass)
    if (!is.na(taxize_lin[r]) && taxize_lin[r] != "Unclassified")
      lineage[r] <- taxize_lin[r]

  lineage
}

# =============================================================================
# STEP 5 — HMMER gene-coverage filter + trimming
# =============================================================================

get_hmm_length <- function(hmm_path) {
  lines <- readLines(hmm_path)
  leng_line <- grep("^LENG", lines, value = TRUE)
  if (length(leng_line) == 0) return(NULL)
  as.integer(strsplit(trimws(leng_line[1]), "\\s+")[[1]][2])
}

write_6frame_fasta <- function(seq_df, out_path) {
  # seq_df: data.frame with columns accession, sequence
  # Writes all 6 reading-frame protein translations to out_path
  # Returns a data.frame: frame_id, accession, strand, offset
  frame_info <- list()
  fh <- file(out_path, "w")
  on.exit(close(fh))

  for (i in seq_len(nrow(seq_df))) {
    acc  <- seq_df$accession[i]
    nseq <- DNAString(seq_df$sequence[i])
    seqs_strands <- list(list(seq = nseq, strand = "+"),
                         list(seq = reverseComplement(nseq), strand = "-"))
    for (ss in seqs_strands) {
      for (offset in 0:2) {
        sub_seq <- subseq(ss$seq, start = offset + 1)
        aa_str  <- as.character(translate(sub_seq, no.init.codon = TRUE,
                                          if.fuzzy.codon = "X"))
        aa_str  <- gsub("\\*", "X", aa_str)
        if (nchar(aa_str) < 10 || nchar(aa_str) > 100000) next

        safe_id  <- gsub("[|;,\\[\\](){}<>\\\\/\\s]", "_", acc, perl = TRUE)
        strand_c <- if (ss$strand == "+") "P" else "M"
        frame_id <- sprintf("%s__F%s%d", safe_id, strand_c, offset)

        cat(sprintf(">%s\n%s\n", frame_id, aa_str), file = fh)
        frame_info[[frame_id]] <- data.frame(
          frame_id  = frame_id,
          accession = acc,
          strand    = ss$strand,
          offset    = offset,
          orig_len  = nchar(seq_df$sequence[i]),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  bind_rows(frame_info)
}

parse_domtblout <- function(domtblout, hmm_len) {
  lines <- readLines(domtblout)
  lines <- lines[!startsWith(lines, "#") & nchar(trimws(lines)) > 0]
  best  <- list()
  for (ln in lines) {
    cols <- strsplit(trimws(ln), "\\s+")[[1]]
    if (length(cols) < 21) next
    frame_id <- cols[1]
    hmm_from <- as.integer(cols[16])
    hmm_to   <- as.integer(cols[17])
    env_from <- as.integer(cols[20])
    env_to   <- as.integer(cols[21])
    coverage <- (hmm_to - hmm_from + 1) / hmm_len
    if (is.null(best[[frame_id]]) || coverage > best[[frame_id]]$coverage) {
      best[[frame_id]] <- list(coverage = coverage, env_from = env_from, env_to = env_to)
    }
  }
  best
}

extract_nt_by_frame <- function(nt_str, seq_len, strand, offset, env_from_aa, env_to_aa) {
  if (strand == "+") {
    nt_start <- max(1, offset + (env_from_aa - 1) * 3 + 1)
    nt_end   <- min(seq_len, offset + env_to_aa * 3)
    return(substr(nt_str, nt_start, nt_end))
  } else {
    rc_start   <- offset + (env_from_aa - 1) * 3
    rc_end     <- offset + env_to_aa * 3
    orig_start <- max(1, seq_len - rc_end + 1)
    orig_end   <- min(seq_len, seq_len - rc_start)
    subseq_str <- substr(nt_str, orig_start, orig_end)
    return(as.character(reverseComplement(DNAString(subseq_str))))
  }
}

run_hmmsearch_filter <- function(seq_df, hmm_profile) {
  hmm_len <- get_hmm_length(hmm_profile)
  if (is.null(hmm_len)) { mcra_log("  Cannot read HMM length", "ERROR"); return(NULL) }

  mcra_log(sprintf("  HMM profile: %s  (length=%d aa, min_cov=%.0f%%)",
                   hmm_profile, hmm_len, HMM_MIN_COVERAGE * 100))

  # Pre-filter by length
  too_long   <- seq_df$accession[nchar(seq_df$sequence) > HMMER_MAX_NT_LEN]
  processable <- seq_df[nchar(seq_df$sequence) <= HMMER_MAX_NT_LEN, ]

  if (length(too_long) > 0)
    mcra_log(sprintf("  %d sequence(s) skipped — length > %d bp", length(too_long), HMMER_MAX_NT_LEN), "WARN")
  if (nrow(processable) == 0) { mcra_log("  No processable sequences", "ERROR"); return(NULL) }

  tmp_dir    <- tempdir()
  tmp_faa    <- file.path(tmp_dir, "mcra_6frames.faa")
  tmp_domtbl <- file.path(tmp_dir, "mcra_hmmsearch.domtbl")

  frame_info <- write_6frame_fasta(processable, tmp_faa)
  mcra_log(sprintf("  %d translated frames written", nrow(frame_info)))

  cmd <- sprintf(
    "%s --domtblout %s --cpu %d --notextw --domE %s %s %s > /dev/null 2>&1",
    HMMSEARCH_PATH, shQuote(tmp_domtbl), HMMER_THREADS, HMMER_EVALUE,
    shQuote(hmm_profile), shQuote(tmp_faa)
  )
  mcra_log(sprintf("  Running: %s", cmd))
  ret <- system(cmd)
  if (ret != 0) { mcra_log("  hmmsearch failed", "ERROR"); return(NULL) }

  raw_hits <- parse_domtblout(tmp_domtbl, hmm_len)
  mcra_log(sprintf("  hmmsearch: %d domain hits", length(raw_hits)))

  # Select best hit per original sequence
  best_per_seq <- list()
  for (frame_id in names(raw_hits)) {
    fi  <- frame_info[frame_info$frame_id == frame_id, ]
    if (nrow(fi) == 0) next
    acc  <- fi$accession[1]
    hit  <- raw_hits[[frame_id]]
    if (is.null(best_per_seq[[acc]]) || hit$coverage > best_per_seq[[acc]]$coverage) {
      best_per_seq[[acc]] <- c(hit, list(strand = fi$strand[1],
                                          offset = fi$offset[1],
                                          orig_len = fi$orig_len[1]))
    }
  }

  # Filter and trim
  accepted <- list()
  n_rejected <- length(too_long)
  for (i in seq_len(nrow(processable))) {
    acc     <- processable$accession[i]
    nt_str  <- processable$sequence[i]
    hit     <- best_per_seq[[acc]]
    if (is.null(hit)) { n_rejected <- n_rejected + 1; next }
    if (hit$coverage < HMM_MIN_COVERAGE) { n_rejected <- n_rejected + 1; next }
    trimmed <- extract_nt_by_frame(nt_str, nchar(nt_str),
                                   hit$strand, hit$offset,
                                   hit$env_from, hit$env_to)
    accepted[[acc]] <- trimmed
  }
  mcra_log(sprintf("  HMMER: %d accepted, %d rejected", length(accepted), n_rejected))
  accepted  # named list: accession → trimmed_seq_str
}


# =============================================================================
# STEP 5b — Template taxonomy annotation transfer
# =============================================================================

load_template_taxonomy <- function(tax_file) {
  # Parse tax4mcrA.taxonomy (tab-separated: accession → lineage string).
  # Returns a named list: accession → named character vector (TAX_LEVELS).
  if (!file.exists(tax_file)) {
    mcra_log(sprintf("  Template taxonomy file not found: '%s' — ",
                     "template annotation transfer skipped", tax_file), "WARN")
    return(list())
  }

  lines <- readLines(tax_file)
  template_tax <- list()

  for (ln in lines) {
    ln <- trimws(ln)
    if (!nzchar(ln)) next
    parts <- strsplit(ln, "\t")[[1]]
    if (length(parts) < 2) next

    accession   <- trimws(parts[1])
    lineage_str <- trimws(gsub(";$", "", parts[2]))   # strip trailing ;
    levels      <- trimws(strsplit(lineage_str, ";")[[1]])

    # Pad to 7 levels
    while (length(levels) < length(TAX_LEVELS)) {
      levels <- c(levels, "Unclassified")
    }
    lineage <- setNames(
      ifelse(nzchar(levels[seq_along(TAX_LEVELS)]),
             levels[seq_along(TAX_LEVELS)],
             "Unclassified"),
      TAX_LEVELS
    )
    template_tax[[accession]] <- lineage
  }

  mcra_log(sprintf("  Loaded %d entries from template taxonomy '%s'",
                   length(template_tax), tax_file))
  template_tax
}


build_blast_db <- function(fasta_path, db_prefix, label = "DB") {
  # Generic makeblastdb wrapper used by both template (Step 5b) and Step 5d.
  # Resolves the binary location, validates inputs, captures stderr, and
  # automatically retries without -parse_seqids if the first attempt fails.

  # ── 1. Resolve makeblastdb binary ───────────────────────────────────────────
  bin <- Sys.which(MAKEBLASTDB_PATH)
  if (!nzchar(bin)) {
    # Try common conda / brew locations as fallback
    candidates <- c(
      file.path(Sys.getenv("CONDA_PREFIX"), "bin", "makeblastdb"),
      "/usr/local/bin/makeblastdb",
      "/usr/bin/makeblastdb"
    )
    bin <- Filter(file.exists, candidates)[1L]
    bin <- if (length(bin) == 1L && !is.na(bin)) bin else ""
  }
  if (!nzchar(bin)) {
    mcra_log(sprintf(
      "  makeblastdb not found (MAKEBLASTDB_PATH='%s'). Set the correct path in config.",
      MAKEBLASTDB_PATH), "WARN")
    return(FALSE)
  }
  mcra_log(sprintf("  makeblastdb binary: %s", bin))

  # ── 2. Validate input FASTA ──────────────────────────────────────────────────
  abs_fasta <- normalizePath(fasta_path, mustWork = FALSE)
  if (!file.exists(abs_fasta)) {
    mcra_log(sprintf("  Input FASTA not found: %s", abs_fasta), "WARN")
    return(FALSE)
  }
  mcra_log(sprintf("  Input FASTA: %s", abs_fasta))

  # ── 3. Attempt makeblastdb (with -parse_seqids first, then without) ──────────
  stderr_f <- tempfile(fileext = ".err")

  run_cmd <- function(extra_flags) {
    cmd <- sprintf(
      "%s -in %s -dbtype nucl -out %s %s 2>%s",
      shQuote(bin), shQuote(abs_fasta), shQuote(db_prefix),
      extra_flags, shQuote(stderr_f)
    )
    mcra_log(sprintf("  Command: %s", cmd))
    system(cmd, ignore.stdout = TRUE)
  }

  ret <- run_cmd("-parse_seqids")
  if (ret != 0) {
    err_txt <- if (file.exists(stderr_f) && file.info(stderr_f)$size > 0L)
                 paste(readLines(stderr_f, warn = FALSE), collapse = " | ")
               else "(no stderr captured)"
    mcra_log(sprintf("  makeblastdb -parse_seqids failed (exit %d): %s",
                     ret, substr(err_txt, 1L, 500L)), "WARN")
    mcra_log("  Retrying without -parse_seqids ...", "WARN")

    ret2 <- run_cmd("")
    if (ret2 != 0) {
      err_txt2 <- if (file.exists(stderr_f) && file.info(stderr_f)$size > 0L)
                    paste(readLines(stderr_f, warn = FALSE), collapse = " | ")
                  else "(no stderr captured)"
      mcra_log(sprintf("  makeblastdb retry failed (exit %d): %s",
                       ret2, substr(err_txt2, 1L, 500L)), "WARN")
      mcra_log(sprintf("  makeblastdb for %s failed — step skipped", label), "WARN")
      return(FALSE)
    }
    mcra_log("  makeblastdb succeeded without -parse_seqids")
  }

  n_seqs <- sum(startsWith(readLines(abs_fasta, warn = FALSE), ">"))
  mcra_log(sprintf("  %s BLAST DB: %s  (%d sequences)", label, db_prefix, n_seqs))
  TRUE
}

build_template_blast_db <- function(template_fasta, db_prefix) {
  build_blast_db(template_fasta, db_prefix, label = "Template")
}


batch_blast_against_template <- function(candidates, template_db_prefix) {
  # candidates: named character vector  accession → seq_str
  if (length(candidates) == 0) return(character(0))

  tmp_dir    <- tempdir()
  query_fa   <- file.path(tmp_dir, "tmpl_query.fa")
  result_tsv <- file.path(tmp_dir, "tmpl_blast.tsv")

  # Safe IDs
  safe_ids <- gsub("[|;,\\[\\](){}<>\\\\/\\s]", "_", names(candidates), perl = TRUE)
  safe_to_orig <- setNames(names(candidates), safe_ids)
  write_fasta(candidates, query_fa, headers = safe_ids)

  cmd <- sprintf(
    paste("%s -query %s -db %s -out %s",
          "-outfmt '6 qseqid sseqid pident length qlen'",
          "-perc_identity %d -max_target_seqs 1 -num_threads %d -dust no"),
    BLASTN_PATH, shQuote(query_fa), shQuote(template_db_prefix),
    shQuote(result_tsv),
    as.integer(TEMPLATE_IDENTITY_THRESHOLD), BLAST_THREADS
  )
  mcra_log(sprintf("  Running template BLAST: %s", cmd))
  ret <- system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
  if (ret != 0 || !file.exists(result_tsv)) {
    mcra_log("  blastn failed for template — template transfer skipped", "WARN")
    return(character(0))
  }

  hits_df <- tryCatch(
    read.table(result_tsv, sep = "\t",
               col.names = c("qseqid","sseqid","pident","length","qlen"),
               stringsAsFactors = FALSE),
    error = function(e) data.frame()
  )
  if (nrow(hits_df) == 0) return(character(0))

  hits_df <- hits_df %>%
    mutate(coverage = length / qlen) %>%
    filter(pident >= TEMPLATE_IDENTITY_THRESHOLD, coverage >= TEMPLATE_QUERY_COVERAGE)

  template_hits <- setNames(hits_df$sseqid, safe_to_orig[hits_df$qseqid])
  template_hits <- template_hits[!is.na(names(template_hits))]

  mcra_log(sprintf("  Template annotation transfer: %d / %d sequences matched",
                   length(template_hits), length(candidates)))
  template_hits  # named character vector: orig_acc → template_accession
}


# =============================================================================
# STEP 5c — Yang et al. annotation transfer
# =============================================================================

load_yang_database <- function(yang_fasta, yang_tax_file) {
  # Accepts the same two-file format as the template step:
  #   yang_fasta    — FASTA with bare accession headers (e.g. >BX950229)
  #   yang_tax_file — tab-separated: ACCESSION\tK;P;C;O;F;G;S;
  missing <- c(
    if (!file.exists(yang_fasta))    sprintf("FASTA '%s'", yang_fasta),
    if (!file.exists(yang_tax_file)) sprintf("taxonomy '%s'", yang_tax_file)
  )
  if (length(missing) > 0) {
    mcra_log(sprintf("  Yang DB file(s) not found: %s — transfer skipped",
                     paste(missing, collapse = ", ")), "WARN")
    return(list(taxonomy = list(), seqs = character(0)))
  }

  # ── Load taxonomy (accession → named lineage vector) ──────────────────────
  tax_lines <- readLines(yang_tax_file)
  yang_tax  <- list()
  for (ln in tax_lines) {
    ln <- trimws(ln)
    if (!nzchar(ln)) next
    parts <- strsplit(ln, "\t")[[1]]
    if (length(parts) < 2) next
    acc         <- trimws(parts[1])
    lineage_str <- trimws(gsub(";$", "", parts[2]))
    levels      <- trimws(strsplit(lineage_str, ";")[[1]])
    while (length(levels) < length(TAX_LEVELS))
      levels <- c(levels, "Unclassified")
    yang_tax[[acc]] <- setNames(
      ifelse(nzchar(levels[seq_along(TAX_LEVELS)]),
             levels[seq_along(TAX_LEVELS)], "Unclassified"),
      TAX_LEVELS
    )
  }

  # ── Load sequences; keep only those with taxonomy ─────────────────────────
  raw      <- read_fasta(yang_fasta)
  taxonomy <- list()
  seqs     <- character(0)

  for (i in seq_along(raw)) {
    acc     <- names(raw)[i]
    seq_str <- gsub("-", "", raw[i])   # strip alignment gaps
    if (is.null(yang_tax[[acc]])) next # no taxonomy entry → skip
    yang_id          <- sprintf("yang_%06d", i)
    taxonomy[[yang_id]] <- yang_tax[[acc]]
    seqs[yang_id]       <- seq_str
  }
  mcra_log(sprintf("  Loaded %d sequences with taxonomy from Yang DB (%d tax entries read)",
                   length(seqs), length(yang_tax)))
  list(taxonomy = taxonomy, seqs = seqs)
}

build_yang_blast_db <- function(yang_seqs, db_prefix) {
  fasta_path <- paste0(db_prefix, ".fasta")
  write_fasta(yang_seqs, fasta_path, headers = names(yang_seqs))
  build_blast_db(fasta_path, db_prefix, label = "Yang")
}

batch_blast_against_yang <- function(candidates, yang_db_prefix) {
  # candidates: named character vector  accession → seq_str
  if (length(candidates) == 0) return(character(0))

  tmp_dir    <- tempdir()
  query_fa   <- file.path(tmp_dir, "yang_query.fa")
  result_tsv <- file.path(tmp_dir, "yang_blast.tsv")

  # Safe IDs
  safe_ids <- gsub("[|;,\\[\\](){}<>\\\\/\\s]", "_", names(candidates), perl = TRUE)
  safe_to_orig <- setNames(names(candidates), safe_ids)
  write_fasta(candidates, query_fa, headers = safe_ids)

  cmd <- sprintf(
    paste("%s -query %s -db %s -out %s",
          "-outfmt '6 qseqid sseqid pident length qlen'",
          "-perc_identity %d -max_target_seqs 1 -num_threads %d -dust no"),
    BLASTN_PATH, shQuote(query_fa), shQuote(yang_db_prefix),
    shQuote(result_tsv),
    as.integer(YANG_IDENTITY_THRESHOLD), BLAST_THREADS
  )
  mcra_log(sprintf("  Running Yang BLAST: %s", cmd))
  ret <- system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
  if (ret != 0 || !file.exists(result_tsv)) {
    mcra_log("  blastn failed — Yang transfer skipped", "WARN")
    return(character(0))
  }

  hits_df <- tryCatch(
    read.table(result_tsv, sep = "\t",
               col.names = c("qseqid","sseqid","pident","length","qlen"),
               stringsAsFactors = FALSE),
    error = function(e) data.frame()
  )
  if (nrow(hits_df) == 0) return(character(0))

  hits_df <- hits_df %>%
    mutate(coverage = length / qlen) %>%
    filter(pident >= YANG_IDENTITY_THRESHOLD, coverage >= YANG_QUERY_COVERAGE)

  yang_hits <- setNames(hits_df$sseqid, safe_to_orig[hits_df$qseqid])
  yang_hits <- yang_hits[!is.na(names(yang_hits))]

  mcra_log(sprintf("  Yang annotation transfer: %d / %d sequences matched",
                   length(yang_hits), length(candidates)))
  yang_hits  # named character vector: orig_acc → yang_id
}


# =============================================================================
# STEP 6 — Filter, deduplicate, build DADA2 entries
# =============================================================================

passes_filters <- function(seq_str) {
  seq_len <- nchar(seq_str)
  if (seq_len < MIN_LENGTH) return(list(ok = FALSE, reason = sprintf("too short (%d bp)", seq_len)))
  if (seq_len > MAX_LENGTH) return(list(ok = FALSE, reason = sprintf("too long (%d bp)", seq_len)))
  n_ambig <- nchar(gsub("[ACGTacgt]", "", seq_str))
  if (seq_len > 0 && n_ambig / seq_len > MAX_AMBIGUOUS_PCT)
    return(list(ok = FALSE, reason = sprintf("too many ambiguous bases (%.1f%%)", n_ambig/seq_len*100)))
  list(ok = TRUE, reason = "")
}

format_dada2_header <- function(lineage, accession = "") {
  # Format: >Kingdom;Phylum;Class;Order;Family;Genus;Species;Accession;
  parts <- c(lineage[TAX_LEVELS], accession)
  paste0(">", paste(parts, collapse = ";"), ";")
}

taxonomy_score <- function(lineage) {
  sum(lineage[TAX_LEVELS] != "Unclassified")
}


# =============================================================================
# STEP 7 — Alignment-based gap trimming (MAFFT)
# =============================================================================

run_mafft <- function(input_fa, output_fa) {
  cmd <- sprintf(
    "%s --auto --thread %d --quiet %s > %s 2>/dev/null",
    MAFFT_PATH, MAFFT_THREADS, shQuote(input_fa), shQuote(output_fa)
  )
  mcra_log(sprintf("  Running MAFFT: %s", cmd))
  ret <- system(cmd)
  if (ret != 0) {
    mcra_log("  MAFFT failed — gap trimming skipped", "WARN")
    return(FALSE)
  }
  TRUE
}

trim_by_gap_profile <- function(seq_id_to_str,
                                 gap_threshold = GAP_THRESHOLD,
                                 min_coverage  = TRIM_MIN_COVERAGE) {
  if (length(seq_id_to_str) == 0) return(list(kept = seq_id_to_str, removed = list()))

  tmp_dir    <- tempdir()
  input_fa   <- file.path(tmp_dir, "mafft_input.fa")
  aligned_fa <- file.path(tmp_dir, "mafft_aligned.fa")

  # Map safe IDs
  safe_ids  <- gsub("[|;,\\[\\](){}<>\\\\/\\s]", "_", names(seq_id_to_str), perl = TRUE)
  # Resolve collisions
  seen_safe <- character(0)
  for (i in seq_along(safe_ids)) {
    base <- safe_ids[i]; k <- 0
    while (safe_ids[i] %in% seen_safe) { k <- k + 1; safe_ids[i] <- sprintf("%s_%d", base, k) }
    seen_safe <- c(seen_safe, safe_ids[i])
  }
  safe_to_orig <- setNames(names(seq_id_to_str), safe_ids)
  names(seq_id_to_str) <- safe_ids   # temp rename

  write_fasta(seq_id_to_str, input_fa)

  ok <- run_mafft(input_fa, aligned_fa)
  if (!ok) {
    names(seq_id_to_str) <- safe_to_orig[safe_ids]   # restore
    return(list(kept = setNames(seq_id_to_str, names(safe_to_orig)),
                removed = list()))
  }

  # Parse alignment
  aln_raw  <- read_fasta(aligned_fa)
  aln_ids  <- names(aln_raw)
  aln_seqs <- toupper(aln_raw)
  aln_len  <- unique(nchar(aln_seqs))
  if (length(aln_len) != 1) {
    mcra_log("  MAFFT alignment has variable sequence lengths — skipping", "WARN")
    return(list(kept = setNames(seq_id_to_str, safe_to_orig[names(seq_id_to_str)]),
                removed = list()))
  }
  n_seqs <- length(aln_ids)

  # Column-wise gap frequencies
  aln_mat   <- do.call(rbind, strsplit(aln_seqs, ""))
  gap_freqs <- colMeans(aln_mat == "-")

  # Block boundaries
  block_start <- which(gap_freqs < gap_threshold)[1]
  block_end   <- tail(which(gap_freqs < gap_threshold), 1)
  if (is.na(block_start) || is.na(block_end)) {
    mcra_log("  No dense block found — keeping all sequences", "WARN")
    return(list(kept = setNames(seq_id_to_str, safe_to_orig[names(seq_id_to_str)]),
                removed = list()))
  }

  block_dense <- gap_freqs[block_start:block_end] < gap_threshold
  n_dense     <- sum(block_dense)
  mcra_log(sprintf("  MAFFT: %d columns, block cols %d–%d, %d dense columns",
                   aln_len, block_start, block_end, n_dense))

  kept    <- list()
  removed <- list()

  for (sid in aln_ids) {
    orig_id      <- safe_to_orig[sid]
    seq_in_block <- substr(aln_seqs[sid], block_start, block_end)
    block_len    <- nchar(seq_in_block)
    chars        <- strsplit(seq_in_block, "")[[1]]

    first_nt <- which(chars != "-")[1]
    last_nt  <- tail(which(chars != "-"), 1)
    if (is.na(first_nt)) {
      removed[[orig_id]] <- "no non-gap characters in block"
      next
    }

    # end_deficit: number of dense columns after last non-gap position
    end_deficit <- if (last_nt < block_len)
      sum(block_dense[(last_nt + 1):block_len])
    else 0L

    # start_deficit: number of dense columns before first non-gap position
    start_deficit <- if (first_nt > 1)
      sum(block_dense[1:(first_nt - 1)])
    else 0L

    if (end_deficit >= min_coverage) {
      removed[[orig_id]] <- sprintf("ends %d dense cols before block-end", end_deficit)
      next
    }
    if (start_deficit >= min_coverage) {
      removed[[orig_id]] <- sprintf("starts %d dense cols after block-start", start_deficit)
      next
    }

    trimmed_seq <- gsub("-", "", seq_in_block)
    if (nchar(trimmed_seq) == 0) {
      removed[[orig_id]] <- "empty after trimming"
      next
    }
    kept[[orig_id]] <- trimmed_seq
  }

  mcra_log(sprintf("  Gap trimming: %d kept, %d removed", length(kept), length(removed)))
  list(kept = kept, removed = removed)
}


# =============================================================================
# STEP 8 — Write outputs
# =============================================================================

write_correspondence_tsv <- function(acc_corr, output_path) {
  df <- bind_rows(lapply(acc_corr, as.data.frame, stringsAsFactors = FALSE))

  # Ensure required columns exist (back-fill defaults if missing)
  if (!"template_annotated" %in% names(df))
    df$template_annotated <- FALSE
  if (!"yang_annotated" %in% names(df))
    df$yang_annotated <- FALSE

  df$annotation_source  <- ifelse(df$template_annotated, "template",
                            ifelse(df$yang_annotated, "yang",
                            ifelse(df$in_final_db, "ncbi", "")))
  df$template_annotated <- ifelse(df$template_annotated, "yes", "no")
  df$yang_annotated     <- ifelse(df$yang_annotated, "yes", "no")
  df$in_final_db        <- ifelse(df$in_final_db, "yes", "no")

  # Reorder columns
  col_order <- c("accession", "dada2_header", "seq_length",
                 "annotation_source", "template_annotated", "yang_annotated",
                 "in_final_db", "status")
  df <- df[, intersect(col_order, names(df))]

  write_tsv(df, output_path)
  n_in    <- sum(df$in_final_db == "yes")
  n_tmpl  <- sum(df$template_annotated == "yes")
  n_yang  <- sum(df$yang_annotated == "yes")
  mcra_log(sprintf("  Correspondence file: %d rows → %s", nrow(df), output_path))
  mcra_log(sprintf("    In final DB          : %d", n_in))
  mcra_log(sprintf("    Template-annotated   : %d", n_tmpl))
  mcra_log(sprintf("    Yang-annotated       : %d", n_yang))
  mcra_log(sprintf("    NCBI-annotated       : %d", n_in - n_tmpl - n_yang))
}


# =============================================================================
# STEP 9 helper functions — Phylogenetic taxonomy curation
# =============================================================================

# ── Resolve a tool binary: PATH first, then common conda/mamba install roots ──
resolve_binary <- function(name) {
  # 0. Expand ~ in case the user supplied an explicit path like ~/anaconda3/bin/X
  expanded <- path.expand(name)
  if (file.exists(expanded) && file.access(expanded, mode = 1L) == 0L)
    return(expanded)

  # 1. Try whatever is in the current PATH
  found <- Sys.which(name)
  if (nzchar(found)) return(found)

  # 2. Try lowercase variant (FastTree vs fasttree)
  lc <- tolower(name)
  if (lc != name) {
    found <- Sys.which(lc)
    if (nzchar(found)) return(found)
  }

  # 3. Search common conda / mamba base-env locations
  home <- Sys.getenv("HOME")
  conda_roots <- c(
    file.path(home, "opt", "anaconda3"),
    file.path(home, "opt", "miniconda3"),
    file.path(home, "anaconda3"),
    file.path(home, "miniconda3"),
    file.path(home, "mambaforge"),
    file.path(home, "miniforge3"),
    file.path(home, "micromamba"),
    "/opt/conda",
    "/opt/anaconda3",
    "/opt/miniconda3"
  )
  for (root in conda_roots) {
    for (nm in unique(c(name, lc))) {
      candidate <- file.path(root, "bin", nm)
      if (file.exists(candidate) && file.access(candidate, mode = 1L) == 0L)
        return(candidate)
      # Also check envs/base/bin for some setups
      candidate2 <- file.path(root, "envs", "base", "bin", nm)
      if (file.exists(candidate2) && file.access(candidate2, mode = 1L) == 0L)
        return(candidate2)
    }
  }

  # 4. Last resort: try `conda run` (slower but guaranteed if conda is on PATH)
  conda_bin <- Sys.which("conda")
  if (nzchar(conda_bin)) {
    for (nm in unique(c(name, lc))) {
      test_cmd <- sprintf("%s run -n base which %s 2>/dev/null",
                          shQuote(conda_bin), nm)
      result <- trimws(system(test_cmd, intern = TRUE, ignore.stderr = TRUE))
      if (length(result) > 0L && nzchar(result[1L]) && file.exists(result[1L]))
        return(result[1L])
    }
  }

  name   # give back the original name and let system() fail with a clear message
}

# ── Build tree: FastTree (GTR) primary, ape NJ fallback ──────────────────────
build_mcra_tree <- function(seq_pool, tree_path) {
  # seq_pool: named character vector  safe_id → ungapped sequence
  # Writes Newick to tree_path.  Returns TRUE on success.

  tmp_dir     <- tempdir()
  input_fa    <- file.path(tmp_dir, "phylo_unaligned.fasta")
  aligned_fa  <- file.path(tmp_dir, "phylo_aligned.fasta")

  write_fasta(seq_pool, input_fa)

  # Align with MAFFT first (required for both FastTree and NJ distance matrix)
  cmd_aln <- sprintf(
    "%s --auto --thread %d --quiet %s > %s 2>/dev/null",
    MAFFT_PATH, MAFFT_THREADS, shQuote(input_fa), shQuote(aligned_fa)
  )
  aln_ok <- (system(cmd_aln) == 0) && file.exists(aligned_fa) &&
            file.info(aligned_fa)$size > 0L
  if (!aln_ok) {
    mcra_log("  MAFFT alignment for tree failed — phylo curation skipped", "WARN")
    return(FALSE)
  }
  mcra_log(sprintf("  Alignment complete (%d sequences)", length(seq_pool)))

  # ── Try FastTree ─────────────────────────────────────────────────────────────
  ft_bin <- resolve_binary(FASTTREE_PATH)
  mcra_log(sprintf("  FastTree binary resolved to: %s", ft_bin))

  cmd_ft <- sprintf(
    "%s -nt -gtr -quiet %s > %s 2>/dev/null",
    shQuote(ft_bin), shQuote(aligned_fa), shQuote(tree_path)
  )
  ft_ok <- (system(cmd_ft) == 0) && file.exists(tree_path) &&
           file.info(tree_path)$size > 0L
  if (ft_ok) {
    mcra_log(sprintf("  FastTree (GTR) complete → %s", tree_path))
    return(TRUE)
  }
  mcra_log(sprintf("  FastTree (%s) not found or failed — trying NJ (ape)", ft_bin), "WARN")

  # ── Fallback: Neighbour-Joining via ape ──────────────────────────────────────
  tryCatch({
    aln_seqs <- read_fasta(aligned_fa)
    mat_list <- lapply(strsplit(toupper(aln_seqs), ""), function(x) x)
    dna_bin  <- ape::as.DNAbin(do.call(rbind, mat_list))
    dm       <- ape::dist.dna(dna_bin, model = "K80", pairwise.deletion = TRUE)
    nj_tree  <- ape::nj(dm)
    ape::write.tree(nj_tree, file = tree_path)
    mcra_log(sprintf("  NJ tree (K80) complete → %s", tree_path))
    TRUE
  }, error = function(e) {
    mcra_log(sprintf("  NJ tree failed: %s", e$message), "WARN")
    FALSE
  })
}

# ── Label each internal node with agreed taxon per rank + distance gate ───────
label_internal_nodes_by_taxonomy <- function(tree, tip_lineages,
                                              patristic_mat, rank_thresholds) {
  n_tips  <- length(tree$tip.label)
  n_nodes <- tree$Nnode
  # node_labels[[i]] = named character vector (TAX_LEVELS) or NULL
  node_labels <- vector("list", n_tips + n_nodes)

  # Pre-compute for each rank: taxon_name → integer vector of tip indices.
  # Used for the monophyletic bifurcation check: a node can be labelled with
  # taxon X at rank R only if ALL tips in the whole tree carrying taxon X at
  # rank R fall within this node's subtree (true monophyly).
  taxon_to_tips <- setNames(vector("list", length(TAX_LEVELS)), TAX_LEVELS)
  for (rk in TAX_LEVELS) {
    tip_map <- list()
    for (i in seq_len(n_tips)) {
      lin <- tip_lineages[[tree$tip.label[i]]]
      if (!is.null(lin) && !is.na(lin[rk]) && lin[rk] != "Unclassified") {
        taxon <- lin[rk]
        tip_map[[taxon]] <- c(tip_map[[taxon]], i)
      }
    }
    taxon_to_tips[[rk]] <- tip_map
  }

  for (node in seq(n_tips + 1L, n_tips + n_nodes)) {
    desc_idx  <- phangorn::Descendants(tree, node, type = "tips")[[1]]
    desc_idx  <- desc_idx[desc_idx <= n_tips]
    if (length(desc_idx) == 0L) next
    tip_names <- tree$tip.label[desc_idx]

    # Maximum patristic distance within this candidate clade
    if (length(tip_names) < 2L) {
      max_pd <- 0
    } else {
      sub_mat <- patristic_mat[tip_names, tip_names, drop = FALSE]
      max_pd  <- max(sub_mat, na.rm = TRUE)
    }

    labels <- setNames(rep(NA_character_, length(TAX_LEVELS)), TAX_LEVELS)

    for (rank in TAX_LEVELS) {
      thr <- rank_thresholds[rank]
      if (is.na(thr) || max_pd > thr) next        # distance gate

      # Collect taxon values from classified tips inside this subtree
      known_taxa <- character(0)
      for (tn in tip_names) {
        lin <- tip_lineages[[tn]]
        if (!is.null(lin) && !is.na(lin[rank]) && lin[rank] != "Unclassified")
          known_taxa <- c(known_taxa, lin[rank])
      }
      if (length(known_taxa) == 0L) next           # no classified members

      uniq <- unique(known_taxa)
      if (length(uniq) != 1L) next                 # polyphyletic inside clade

      taxon_X <- uniq[1L]

      # ── True monophyletic bifurcation check ──────────────────────────────────
      # ALL tips in the full tree carrying taxon_X at this rank must be
      # contained within this subtree.  If any lie outside, the group is
      # not monophyletic and this node must not be labelled.
      all_tips_with_X <- taxon_to_tips[[rank]][[taxon_X]]
      if (!all(all_tips_with_X %in% desc_idx)) next   # paraphyletic in full tree

      labels[rank] <- taxon_X
    }

    if (any(!is.na(labels))) node_labels[[node]] <- labels
  }

  node_labels
}

# ── Walk each unclassified tip root-ward; assign deepest agreeing label ───────
curate_taxonomy_by_phylo <- function(tree, tip_lineages, node_labels) {
  n_tips    <- length(tree$tip.label)
  curated   <- tip_lineages
  n_updated <- 0L
  detail    <- list()   # safe_id → list(before, after)

  for (tip_idx in seq_len(n_tips)) {
    tn  <- tree$tip.label[tip_idx]
    lin <- tip_lineages[[tn]]
    if (is.null(lin)) next

    unknown_ranks <- TAX_LEVELS[!is.na(lin[TAX_LEVELS]) & lin[TAX_LEVELS] == "Unclassified"]
    if (length(unknown_ranks) == 0L) next   # fully classified

    # phangorn::Ancestors returns nodes in tip-to-root order:
    #   ancestors[1] = immediate parent (closest to tip, most specific)
    #   ancestors[k] = k-th ancestor (further from tip, less specific)
    ancestors <- phangorn::Ancestors(tree, tip_idx, type = "all")
    if (length(ancestors) == 0L) next

    new_lin     <- lin
    changed     <- FALSE
    # max_anc_idx: the highest ancestor index (1 = closest to tip) eligible
    # for assigning the current and all subsequent ranks.  Starts at the full
    # ancestor chain; shrinks each time a rank is assigned, ensuring that
    # finer ranks (e.g. Family) can only come from nodes that are between the
    # tip and the node that provided the coarser rank (e.g. Order) — top-down
    # consistency.  Without this constraint, different ranks could be assigned
    # from incompatible branches of the tree.
    max_anc_idx <- length(ancestors)

    for (rank in TAX_LEVELS) {          # TAX_LEVELS = coarse → fine order
      if (!rank %in% unknown_ranks) next # already annotated; skip
      if (max_anc_idx < 1L) break        # no eligible ancestors remain

      for (k in seq_len(max_anc_idx)) {
        anc_lbl <- node_labels[[ancestors[k]]]
        if (!is.null(anc_lbl) && !is.na(anc_lbl[rank])) {
          new_lin[rank] <- anc_lbl[rank]
          changed       <- TRUE
          max_anc_idx   <- k - 1L   # finer ranks restricted to [1, k-1]
          break
        }
      }
      # If no label found for this rank within the eligible window, we do NOT
      # break — a finer rank may still be assignable from within the same window.
    }

    if (changed) {
      curated[[tn]] <- new_lin
      n_updated     <- n_updated + 1L
      detail[[tn]]  <- list(before = lin, after = new_lin)
    }
  }

  list(lineages = curated, n_updated = n_updated, detail = detail)
}

# ── Leaf-to-root curation with consensus lineage fill ─────────────────────────
#
# For each unclassified tip, walk from the tip toward the root.
# At the FIRST (deepest = closest to tip) ancestor that carries a monophyletic
# label within the patristic distance threshold:
#   1. Assign that rank label to the tip.
#   2. For every COARSER rank that is still Unclassified, derive the value by
#      majority consensus of the known sequences in that ancestor's subtree.
#      This fills the complete lineage ("Methanococcaceae" → Order, Class,
#      Phylum, Domain are filled from the subtree's known annotations).
#
# Unlike the top-down variant, the anchor is always the deepest available node
# so the most specific monophyletic clade is preferred.
curate_taxonomy_leaf_to_root <- function(tree, tip_lineages, node_labels) {
  n_tips    <- length(tree$tip.label)
  curated   <- tip_lineages
  n_updated <- 0L
  detail    <- list()

  # Pre-compute subtree tip lists for every internal node (reused per tip walk)
  subtree_tips <- lapply(seq(n_tips + 1L, n_tips + tree$Nnode), function(nd) {
    idx <- phangorn::Descendants(tree, nd, type = "tips")[[1L]]
    tree$tip.label[idx[idx <= n_tips]]
  })
  names(subtree_tips) <- as.character(seq(n_tips + 1L, n_tips + tree$Nnode))

  # Consensus value for a rank across a set of known tips
  rank_consensus <- function(tips, rank) {
    vals <- vapply(tips, function(t) {
      l <- tip_lineages[[t]]
      if (!is.null(l) && !is.na(l[rank]) && l[rank] != "Unclassified") l[rank]
      else NA_character_
    }, character(1L))
    vals <- vals[!is.na(vals)]
    uniq <- unique(vals)
    if (length(uniq) == 1L) uniq else NA_character_
  }

  for (tip_idx in seq_len(n_tips)) {
    tn  <- tree$tip.label[tip_idx]
    lin <- tip_lineages[[tn]]
    if (is.null(lin)) next

    unknown_ranks <- TAX_LEVELS[!is.na(lin[TAX_LEVELS]) & lin[TAX_LEVELS] == "Unclassified"]
    if (length(unknown_ranks) == 0L) next

    ancestors <- phangorn::Ancestors(tree, tip_idx, type = "all")
    if (length(ancestors) == 0L) next

    new_lin <- lin
    changed <- FALSE

    # Walk tip → root; stop at the FIRST ancestor with any rank label
    for (k in seq_along(ancestors)) {
      anc     <- ancestors[k]
      anc_lbl <- node_labels[[anc]]
      if (is.null(anc_lbl) || all(is.na(anc_lbl))) next

      # Assign all labeled ranks from this anchor node
      for (rank in unknown_ranks) {
        if (!is.na(anc_lbl[rank])) {
          new_lin[rank] <- anc_lbl[rank]
          changed <- TRUE
        }
      }

      # For still-Unclassified ranks (typically coarser ranks above the anchor),
      # derive by consensus of known sequences in this ancestor's subtree.
      subtree_key  <- as.character(anc)
      anc_tips     <- subtree_tips[[subtree_key]]
      if (!is.null(anc_tips) && length(anc_tips) > 0L) {
        for (rank in unknown_ranks) {
          if (new_lin[rank] != "Unclassified") next   # already assigned
          consensus_val <- rank_consensus(anc_tips, rank)
          if (!is.na(consensus_val)) {
            new_lin[rank] <- consensus_val
            changed <- TRUE
          }
        }
      }

      break   # only use the deepest anchor found
    }

    if (changed) {
      curated[[tn]] <- new_lin
      n_updated     <- n_updated + 1L
      detail[[tn]]  <- list(before = lin, after = new_lin)
    }
  }

  list(lineages = curated, n_updated = n_updated, detail = detail)
}


# =============================================================================
# STEP 9a / 9c helper functions — neighbour-guided label cleaning and filling
# =============================================================================

UNK <- "Unclassified"

changes <- list()
record_changes <- function(stage, i, k, old, new, n_top = NA, n_total = NA) {
  if (!length(i)) return(invisible())
  changes[[length(changes) + 1L]] <<- data.frame(
    stage = stage, accession = db$acc[i], rank = RANKS[k], old = old, new = new,
    n_support = n_top, n_voting = n_total, stringsAsFactors = FALSE)
}

unclassified_counts <- function(lab, stage) {
  data.frame(stage = stage, rank = RANKS,
             unclassified = colSums(lab[, RANKS, drop = FALSE] == UNK),
             distinct_names = vapply(RANKS, function(r) length(setdiff(unique(lab[, r]), UNK)), 0L),
             row.names = NULL)
}

# All-vs-all nearest neighbours (self hits excluded by unique labels + --self)
run_vsearch_neighbours <- function(seqs, out_tsv) {
  fa <- tempfile(fileext = ".fasta")
  writeLines(paste0(">", seq_along(seqs), "\n", seqs), fa)
  cmd <- sprintf(paste("%s --usearch_global %s --db %s --self --strand plus",
                       "--id %.2f --maxaccepts %d --maxrejects %d --threads %d",
                       "--userfields query+target+id --userout %s --quiet"),
                 VSEARCH_PATH, shQuote(fa), shQuote(fa), NEIGHBOUR_MIN_ID,
                 NEIGHBOUR_MAX_HITS, NEIGHBOUR_MAX_HITS * 10L, VSEARCH_THREADS, shQuote(out_tsv))
  if (system(cmd) != 0) stop("vsearch failed: ", cmd)
  h <- read.delim(out_tsv, header = FALSE, col.names = c("q", "t", "id"))
  h$id <- h$id / 100
  h
}

# Youden J cut-off on % identity: same-taxon vs different-taxon neighbour pairs
calibrate_cutoffs <- function(hits, lab, voter) {
  grid <- seq(NEIGHBOUR_MIN_ID, 1, by = 0.001)
  res <- lapply(seq_along(RANKS), function(k) {
    ok <- voter[hits$q, k] & voter[hits$t, k] &
          lab[hits$q, k] != UNK & lab[hits$t, k] != UNK
    same <- lab[hits$q[ok], k] == lab[hits$t[ok], k]
    id   <- hits$id[ok]
    n_pos <- sum(same); n_neg <- sum(!same)
    if (n_pos < CALIB_MIN_PAIRS || n_neg < CALIB_MIN_PAIRS)
      return(data.frame(rank = RANKS[k], youden_cutoff = NA, J = NA, precision_cutoff = NA,
                        n_same = n_pos, n_diff = n_neg))
    pos <- sort(id[same]); neg <- sort(id[!same])
    sens <- 1 - findInterval(grid, pos, left.open = TRUE) / n_pos   # P(id >= t | same)
    spec <- findInterval(grid, neg, left.open = TRUE) / n_neg       # P(id <  t | diff)
    J <- sens + spec - 1
    b <- which.max(J)
    # Precision cut-off: lowest identity t such that, for every s >= t, at
    # least PRECISION_TARGET of the pairs with identity >= s share the taxon
    n_above_pos <- n_pos - findInterval(grid, pos, left.open = TRUE)
    n_above_neg <- n_neg - findInterval(grid, neg, left.open = TRUE)
    prec <- n_above_pos / pmax(1, n_above_pos + n_above_neg)
    okp  <- rev(cumprod(rev(prec >= PRECISION_TARGET))) == 1 &
            (n_above_pos + n_above_neg) >= CALIB_MIN_PAIRS
    pc <- if (any(okp)) grid[which(okp)[1L]] else NA
    data.frame(rank = RANKS[k], youden_cutoff = grid[b], J = round(J[b], 3),
               precision_cutoff = pc, n_same = n_pos, n_diff = n_neg)
  })
  res <- do.call(rbind, res)
  # The precision cut-off is used (Youden J is reported for comparison).
  # A rank without enough pairs takes the cut-off of the nearest finer
  # calibrated rank (sharing a finer taxon implies sharing the coarser one).
  # Cut-offs are then made non-decreasing from kingdom to genus, so a finer
  # rank never accepts neighbours that a coarser rank would reject.
  cut <- res$precision_cutoff
  for (k in rev(seq_along(cut))) if (is.na(cut[k]))
    cut[k] <- if (k < length(cut)) cut[k + 1L] else res$youden_cutoff[k]
  if (anyNA(cut)) stop("no identity cut-off could be calibrated for: ",
                       paste(RANKS[is.na(cut)], collapse = ", "))
  res$cutoff <- cummax(cut)
  res$source <- ifelse(is.na(res$precision_cutoff), "nearest finer rank",
                ifelse(res$cutoff > res$precision_cutoff, "raised to coarser rank", "precision"))
  # Share of voter pairs at or above the cut-off that carry the same taxon
  res$same_taxon_above_cutoff <- vapply(seq_along(RANKS), function(k) {
    ok <- voter[hits$q, k] & voter[hits$t, k] & lab[hits$q, k] != UNK &
          lab[hits$t, k] != UNK & hits$id >= res$cutoff[k]
    round(mean(lab[hits$q[ok], k] == lab[hits$t[ok], k]), 4)
  }, 0)
  res
}

# Majority vote of usable neighbours at rank k for every query
vote_rank <- function(k, lab, voter, hits, cutoff) {
  h <- hits[hits$id >= cutoff, c("q", "t")]
  v <- lab[h$t, k]
  ok <- voter[h$t, k] & v != UNK
  h <- h[ok, ]; v <- v[ok]
  if (!nrow(h)) return(NULL)
  key <- paste(h$q, v, sep = "\t")
  tab <- table(key)
  sp  <- strsplit(names(tab), "\t", fixed = TRUE)
  df  <- data.frame(q = as.integer(vapply(sp, `[`, "", 1L)),
                    taxon = vapply(sp, `[`, "", 2L), n = as.integer(tab),
                    stringsAsFactors = FALSE)
  df$total <- ave(df$n, df$q, FUN = sum)
  df  <- df[order(df$q, -df$n), ]
  top <- df[!duplicated(df$q), ]
  own_key <- paste(top$q, lab[top$q, k], sep = "\t")
  top$n_own <- as.integer(tab[own_key]); top$n_own[is.na(top$n_own)] <- 0L
  top
}

# name at rank k → parent name at rank k-1 supported by >= MIN_AGREEMENT of voters
canonical_parents <- function(lab, voter) {
  out <- vector("list", length(RANKS))
  for (k in 2:length(RANKS)) {
    ok <- voter[, k] & voter[, k - 1L] & lab[, k] != UNK & lab[, k - 1L] != UNK
    if (!any(ok)) { out[[k]] <- character(0); next }
    tb <- table(lab[ok, k], lab[ok, k - 1L])
    best <- apply(tb, 1L, function(x) if (max(x) / sum(x) >= MIN_AGREEMENT)
                  colnames(tb)[which.max(x)] else NA_character_)
    out[[k]] <- best[!is.na(best)]
  }
  out
}

canonical_lineage <- function(name, k, parents) {
  lin <- rep(NA_character_, k); lin[k] <- name
  if (k > 1L) for (j in k:2) {
    p <- parents[[j]][lin[j]]
    if (is.na(p)) break
    lin[j - 1L] <- p
  }
  lin
}

# Each named rank takes its canonical parent (fills "Unclassified" parents too)
enforce_lineages <- function(lab, parents, protect, stage) {
  for (k in length(RANKS):2) {
    for (i in which(lab[, k] != UNK & !protect)) {
      p <- parents[[k]][lab[i, k]]
      if (is.na(p) || lab[i, k - 1L] == p) next
      record_changes(stage, i, k - 1L, lab[i, k - 1L], p)
      lab[i, k - 1L] <- p
    }
  }
  lab
}

# A4 / C1: replace or blank named labels contradicted by near neighbours
check_named_labels <- function(lab, src, voter, hits, cut, parents,
                               candidates, stage, new_src) {
  for (k in seq_along(RANKS)) {
    top <- vote_rank(k, lab, voter, hits, cut[k])
    if (is.null(top)) next
    top <- top[candidates[top$q, k] & lab[top$q, k] != UNK, ]
    fix <- top[top$taxon != lab[top$q, k] & top$n >= MIN_SUPPORT_FIX &
               top$n / top$total >= MIN_AGREEMENT, ]
    blank <- top[top$n_own == 0L & top$total >= MIN_SUPPORT_FIX &
                 !(top$q %in% fix$q), ]
    for (r in seq_len(nrow(fix))) {
      i <- fix$q[r]; Y <- fix$taxon[r]
      record_changes(stage, i, k, lab[i, k], Y, fix$n[r], fix$total[r])
      lab[i, k] <- Y; src[i, k] <- new_src
      canon <- canonical_lineage(Y, k, parents)
      for (j in seq_len(k - 1L)) if (!is.na(canon[j]) && lab[i, j] != canon[j]) {
        record_changes(stage, i, j, lab[i, j], canon[j]); lab[i, j] <- canon[j]; src[i, j] <- new_src
      }
      if (k < length(RANKS)) for (j in (k + 1L):length(RANKS)) {
        if (lab[i, j] == UNK) next
        if (!identical(unname(canonical_lineage(lab[i, j], j, parents)[k]), Y)) {
          record_changes(stage, i, j, lab[i, j], UNK)
          lab[i, j:length(RANKS)] <- UNK; src[i, j:length(RANKS)] <- new_src
          break
        }
      }
    }
    for (r in seq_len(nrow(blank))) {
      i <- blank$q[r]
      record_changes(paste0(stage, "_blank"), i, k, lab[i, k], UNK, blank$n[r], blank$total[r])
      lab[i, k:length(RANKS)] <- UNK; src[i, k:length(RANKS)] <- new_src
    }
    mcra_log(sprintf("  %-8s cut-off %.3f: %d replaced, %d blanked",
                     RANKS[k], cut[k], nrow(fix), nrow(blank)))
  }
  list(lab = lab, src = src)
}

# C2: fill "Unclassified" ranks from agreeing original-record neighbours
fill_unclassified <- function(lab, src, voter, hits, cut, parents, blind, stage) {
  for (k in seq_along(RANKS)) {
    top <- vote_rank(k, lab, voter, hits, cut[k])
    if (is.null(top)) next
    top <- top[lab[top$q, k] == UNK & !blind[top$q] &
               top$n >= MIN_SUPPORT_FILL & top$n / top$total >= MIN_AGREEMENT, ]
    if (k > 1L) {
      par <- parents[[k]][top$taxon]
      top <- top[!is.na(par) & lab[top$q, k - 1L] == par, ]
    }
    record_changes(stage, top$q, k, rep(UNK, nrow(top)), top$taxon, top$n, top$total)
    lab[cbind(top$q, k)] <- top$taxon
    src[cbind(top$q, k)] <- "neighbour_fill"
    mcra_log(sprintf("  %-8s cut-off %.3f: %d filled", RANKS[k], cut[k], nrow(top)))
  }
  list(lab = lab, src = src)
}

# =============================================================================
# STEP 9 node labelling — exact, memory-light re-implementation
# =============================================================================
# label_internal_nodes_by_taxonomy() builds the full patristic matrix
# (27,942² doubles ≈ 6 GB) and scans every sub-matrix.  This version computes
# the same three conditions in one post-order pass:
#   max patristic distance within a clade = clade diameter
#     (max over internal nodes of the two deepest child paths),
#   unanimity = the clade's classified tips carry one taxon,
#   strict monophyly = the clade holds every tip of that taxon in the tree.
label_internal_nodes_fast <- function(tree, tip_lineages, rank_thresholds) {
  n_tips <- length(tree$tip.label); n_all <- n_tips + tree$Nnode
  tr <- reorder(tree, "postorder")
  bl <- if (is.null(tr$edge.length)) rep(0, nrow(tr$edge)) else tr$edge.length
  height <- numeric(n_all); diam <- numeric(n_all)
  top1 <- rep(-Inf, n_all); top2 <- rep(-Inf, n_all)
  root <- n_tips + 1L
  ranks <- names(rank_thresholds)
  tx  <- matrix(0L, n_all, length(ranks), dimnames = list(NULL, ranks))
  cnt <- matrix(0L, n_all, length(ranks))
  taxa <- list(); total <- list()
  for (r in seq_along(ranks)) {
    v <- vapply(tree$tip.label, function(t) {
      l <- tip_lineages[[t]]
      if (is.null(l) || is.na(l[ranks[r]]) || l[ranks[r]] == "Unclassified") NA_character_
      else unname(l[ranks[r]])
    }, "")
    f <- factor(v); taxa[[r]] <- levels(f)
    tx[seq_len(n_tips), r] <- ifelse(is.na(f), 0L, as.integer(f))
    cnt[seq_len(n_tips), r] <- as.integer(!is.na(f))
    total[[r]] <- tabulate(as.integer(f[!is.na(f)]), nbins = length(levels(f)))
  }
  finalize <- function(p) {
    if (p > n_tips) {
      height[p] <<- top1[p]
      if (is.finite(top2[p])) diam[p] <<- max(diam[p], top1[p] + top2[p])
    }
  }
  for (e in seq_len(nrow(tr$edge))) {
    p <- tr$edge[e, 1L]; c <- tr$edge[e, 2L]
    finalize(c)
    val <- height[c] + bl[e]
    if (val > top1[p]) { top2[p] <- top1[p]; top1[p] <- val }
    else if (val > top2[p]) top2[p] <- val
    diam[p] <- max(diam[p], diam[c])
    for (r in seq_along(ranks)) {
      tc <- tx[c, r]
      if (tc == 0L) next
      tp <- tx[p, r]
      if (tp == 0L) { tx[p, r] <- tc; cnt[p, r] <- cnt[c, r] }
      else if (tp > 0L && tp == tc) cnt[p, r] <- cnt[p, r] + cnt[c, r]
      else tx[p, r] <- -1L
    }
  }
  finalize(root)
  node_labels <- vector("list", n_all)
  for (node in seq(n_tips + 1L, n_all)) {
    labels <- setNames(rep(NA_character_, length(TAX_LEVELS)), TAX_LEVELS)
    for (r in seq_along(ranks)) {
      thr <- rank_thresholds[r]
      if (is.na(thr) || diam[node] > thr) next
      t <- tx[node, r]
      if (t <= 0L) next
      if (cnt[node, r] != total[[r]][t]) next
      labels[ranks[r]] <- taxa[[r]][t]
    }
    if (any(!is.na(labels))) node_labels[[node]] <- labels
  }
  node_labels
}

# curate_taxonomy_leaf_to_root() pre-computes the tip list of every internal
# node (≈ n²/2 entries on a deep 27,942-tip tree).  Only labelled nodes can be
# anchors, so this version lists tips for those nodes alone: each tip's anchor
# is the deepest labelled node above it, and the consensus of each anchor's
# subtree is computed once.  Assignments are the same as the original.
curate_leaf_to_root_fast <- function(tree, tip_lineages, node_labels) {
  n_tips <- length(tree$tip.label)
  curated <- tip_lineages; detail <- list()
  lab_nodes <- which(!vapply(node_labels, is.null, TRUE))
  lab_nodes <- lab_nodes[lab_nodes > n_tips]
  if (!length(lab_nodes)) return(list(lineages = curated, n_updated = 0L, detail = detail))
  n_anc <- vapply(lab_nodes, function(nd) length(phangorn::Ancestors(tree, nd, type = "all")), 0L)
  lab_nodes <- lab_nodes[order(n_anc)]          # root-most first
  anchor <- integer(n_tips)
  sub_tips <- list()
  for (nd in lab_nodes) {
    idx <- phangorn::Descendants(tree, nd, type = "tips")[[1L]]
    idx <- idx[idx <= n_tips]
    sub_tips[[as.character(nd)]] <- idx
    anchor[idx] <- nd                           # deeper labelled nodes overwrite
  }
  lin_mat <- do.call(rbind, lapply(tree$tip.label, function(t) {
    l <- tip_lineages[[t]]; if (is.null(l)) rep(NA_character_, length(TAX_LEVELS)) else l[TAX_LEVELS]
  }))
  colnames(lin_mat) <- TAX_LEVELS
  consensus <- lapply(sub_tips, function(idx) vapply(TAX_LEVELS, function(r) {
    v <- lin_mat[idx, r]; v <- unique(v[!is.na(v) & v != "Unclassified"])
    if (length(v) == 1L) v else NA_character_
  }, ""))
  for (i in which(anchor > 0L)) {
    tn <- tree$tip.label[i]; lin <- tip_lineages[[tn]]
    if (is.null(lin)) next
    unknown <- TAX_LEVELS[!is.na(lin[TAX_LEVELS]) & lin[TAX_LEVELS] == "Unclassified"]
    if (!length(unknown)) next
    a <- as.character(anchor[i]); lbl <- node_labels[[anchor[i]]]; cons <- consensus[[a]]
    new_lin <- lin; changed <- FALSE
    for (r in unknown) if (!is.na(lbl[r])) { new_lin[r] <- lbl[r]; changed <- TRUE }
    for (r in unknown) {
      if (new_lin[r] != "Unclassified") next
      if (!is.na(cons[r])) { new_lin[r] <- cons[r]; changed <- TRUE }
    }
    if (changed) { curated[[tn]] <- new_lin; detail[[tn]] <- list(before = lin, after = new_lin) }
  }
  list(lineages = curated, n_updated = length(detail), detail = detail)
}

run_step9 <- function(lab, tree, blind) {
  safe <- tree$tip.label
  lins <- lapply(seq_len(nrow(lab)), function(i) {
    l <- lab[i, ]
    if (blind[i]) l[RANKS] <- UNK
    l
  })
  names(lins) <- db$safe_id
  node_labels <- if (FAST_NODE_LABELLING)
    label_internal_nodes_fast(tree, lins, PHYLO_RANK_THRESHOLDS)
  else
    label_internal_nodes_by_taxonomy(tree, lins,
      as.matrix(ape::cophenetic.phylo(tree)), PHYLO_RANK_THRESHOLDS)
  mcra_log(sprintf("  %d / %d internal nodes received at least one rank label",
                   sum(!vapply(node_labels, is.null, TRUE)), tree$Nnode))
  res <- if (FAST_NODE_LABELLING) curate_leaf_to_root_fast(tree, lins, node_labels)
         else curate_taxonomy_leaf_to_root(tree, lins, node_labels)
  mcra_log(sprintf("  Leaf-to-root curation: %d tip(s) received new rank labels", res$n_updated))
  if (CORRECT_TAXONOMY && res$n_updated > 0L) {
    mcra_log("  Verifying newly assigned names against NCBI (taxize)...")
    for (sid in names(res$detail))
      res$lineages[[sid]] <- correct_lineage_names(res$lineages[[sid]])
  }
  if (CORRECT_TAXONOMY) {
    mcra_log("  Filling Unclassified higher ranks via taxize lineage lookup...")
    for (sid in names(res$lineages))
      res$lineages[[sid]] <- fill_lineage_from_taxize(res$lineages[[sid]])
  }
  out <- do.call(rbind, lapply(db$safe_id, function(s) res$lineages[[s]][TAX_LEVELS]))
  out[blind, RANKS] <- lab[blind, RANKS]   # blinded tips keep their own labels
  out
}

# =============================================================================
# CURATION PIPELINE (blind = references hidden from every decision)
# =============================================================================

curate <- function(blind = rep(FALSE, length(db$acc)), quiet_stage = "") {
  lab <- db$lab
  # A1 — one spelling for "unknown"
  for (r in RANKS) {
    unres <- grepl(REFINE_UNRES_PATTERN, lab[, r])
    if (!any(blind)) record_changes("A1_placeholder", which(unres & lab[, r] != UNK), match(r, RANKS),
                                    lab[unres & lab[, r] != UNK, r], UNK)
    lab[unres, r] <- UNK
  }
  if (!any(blind)) summary_rows[[length(summary_rows) + 1L]] <<- unclassified_counts(lab, "A1 placeholders collapsed")

  genome <- grepl("^GC[AF]_", db$acc)
  src <- matrix(ifelse(genome, "genome", "original"), length(db$acc), length(RANKS))
  src[lab[, RANKS] == UNK] <- "none"
  voter_of <- function(src) (src %in% c("genome", "original", "neighbour_fix")) & !blind
  voter <- matrix(voter_of(src), ncol = length(RANKS))

  # A2 — identity cut-offs
  cal <- calibrate_cutoffs(hits, lab, voter)
  cut <- setNames(cal$cutoff, RANKS)
  if (!any(blind)) {
    write.table(cal, OUTPUT_FILE_CURATION_CUTOFFS, sep = "\t", quote = FALSE, row.names = FALSE)
    for (r in seq_len(nrow(cal)))
      mcra_log(sprintf("  %-8s identity cut-off %.3f  (%s; precision %s, Youden %s, J=%s, same=%d, diff=%d; same taxon above cut-off %.2f%%)",
                       cal$rank[r], cal$cutoff[r], cal$source[r], cal$precision_cutoff[r], cal$youden_cutoff[r], cal$J[r],
                       cal$n_same[r], cal$n_diff[r], 100 * cal$same_taxon_above_cutoff[r]))
  }

  # A3 — lineage consistency
  parents <- canonical_parents(lab, voter)
  lab <- enforce_lineages(lab, parents, protect = genome | blind, stage = "A3_lineage")

  # A4 — neighbour check of named labels
  mcra_log("[A4] Neighbour check of named labels")
  cand <- matrix(!(genome | blind), length(db$acc), length(RANKS))
  r <- check_named_labels(lab, src, voter, hits, cut, parents, cand, "A4_neighbour_fix", "neighbour_fix")
  lab <- r$lab; src <- r$src
  if (!any(blind)) summary_rows[[length(summary_rows) + 1L]] <<- unclassified_counts(lab, "A pre-Step 9 cleaning")

  # B — Step 9 re-run
  mcra_log("[B] Step 9 re-run")
  lab9 <- run_step9(lab, tree, blind)
  added <- lab[, RANKS] == UNK & lab9[, RANKS] != UNK
  changed <- lab[, RANKS] != UNK & lab9[, RANKS] != lab[, RANKS]
  if (!any(blind)) for (k in seq_along(RANKS)) {
    record_changes("B_step9", which(added[, k]), k, lab[added[, k], k], lab9[added[, k], k])
    record_changes("B_step9_taxize", which(changed[, k]), k, lab[changed[, k], k], lab9[changed[, k], k])
  }
  src[added | changed] <- "step9"
  lab <- lab9
  if (!any(blind)) summary_rows[[length(summary_rows) + 1L]] <<- unclassified_counts(lab, "B Step 9 re-run")

  # C1 — neighbour check of the labels Step 9 added
  mcra_log("[C1] Neighbour check of Step 9 labels")
  voter <- matrix(voter_of(src), ncol = length(RANKS))
  r <- check_named_labels(lab, src, voter, hits, cut, parents, src == "step9" & !blind,
                          "C1_step9_fix", "step9_fix")
  lab <- r$lab; src <- r$src

  # C2 — neighbour fill
  mcra_log("[C2] Neighbour fill of Unclassified ranks")
  r <- fill_unclassified(lab, src, voter, hits, cut, parents, blind, "C2_neighbour_fill")
  lab <- r$lab; src <- r$src

  # C3 — lineage consistency
  lab <- enforce_lineages(lab, parents, protect = genome | blind, stage = "C3_lineage")
  if (!any(blind)) summary_rows[[length(summary_rows) + 1L]] <<- unclassified_counts(lab, "C post-Step 9")
  lab
}

# =============================================================================
# OUTPUT WRITERS
# =============================================================================

write_full_db <- function(lab, path) {
  hdr <- vapply(seq_along(db$acc), function(i) format_dada2_header(lab[i, ], db$acc[i]), "")
  writeLines(paste0(hdr, "\n", db$seq), path)
}

# Kingdom…Genus, truncated at the first unknown rank; exact duplicates collapse
# to their longest common lineage prefix; unknown at TRAIN_MIN_RANK → dropped
training_set <- function(lab, keep = rep(TRUE, nrow(lab))) {
  lab <- lab[keep, RANKS, drop = FALSE]; seqs <- db$seq[keep]
  depth <- apply(lab != UNK, 1L, function(x) if (all(x)) length(x) else which.min(x) - 1L)
  grp <- split(seq_along(seqs), seqs)
  out_seq <- character(0); out_lin <- character(0)
  for (g in grp) {
    d <- min(depth[g])
    if (length(g) > 1L && d > 0L) {
      same <- vapply(seq_len(d), function(k) length(unique(lab[g, k])) == 1L, TRUE)
      d <- if (all(same)) d else which.min(same) - 1L
    }
    if (d < match(TRAIN_MIN_RANK, RANKS)) next
    out_seq <- c(out_seq, seqs[g[1L]])
    out_lin <- c(out_lin, paste0(paste(lab[g[1L], seq_len(d)], collapse = ";"), ";"))
  }
  list(seq = out_seq, lin = out_lin)
}

species_set <- function(lab, keep = rep(TRUE, nrow(lab))) {
  sp <- db$lab[, "species"]
  words <- strsplit(sp, " ", fixed = TRUE)
  g <- vapply(words, `[`, "", 1L); e <- vapply(words, function(w) if (length(w) > 1L) w[2L] else "", "")
  ok <- keep & lab[, "genus"] != UNK & g == lab[, "genus"] &
        grepl("^[a-z][a-z-]+$", e) & !e %in% c("sp", "spp", "archaeon", "cf", "aff")
  list(id = db$acc[ok], seq = db$seq[ok], name = paste(g[ok], e[ok]))
}

# ── Step 9 driver: 9a cleaning → 9b tree curation → 9c fill; writes 9d files ──
refine_and_curate_step9 <- function(accs, safe_pool, safe_lins, phy_tree) {
  RANKS <<- setdiff(TAX_LEVELS, "species")
  lab <- do.call(rbind, lapply(names(safe_pool), function(s) {
    l <- safe_lins[[s]][TAX_LEVELS]; l[is.na(l)] <- UNK; l }))
  colnames(lab) <- TAX_LEVELS
  db <<- list(acc = unname(accs), seq = toupper(unname(safe_pool)), lab = lab,
              safe_id = names(safe_pool))
  tree <<- phy_tree
  changes <<- list()
  summary_rows <<- list(unclassified_counts(lab, "input (Step 8 lineages)"))
  mcra_log("  [9a] Nearest-neighbour search (vsearch, all vs all)")
  hits <<- run_vsearch_neighbours(db$seq, NEIGHBOURS_FILE)
  mcra_log(sprintf("  %d neighbour pairs at >= %.2f identity", nrow(hits), NEIGHBOUR_MIN_ID))
  lab_final <- curate()

  ts <- training_set(lab_final)
  writeLines(paste0(">", ts$lin, "\n", ts$seq), OUTPUT_FILE_CURATED_TRAIN)
  ss <- species_set(lab_final)
  writeLines(paste0(">", ss$id, " ", ss$name, "\n", ss$seq), OUTPUT_FILE_CURATED_SPECIES)
  ch <- if (length(changes)) do.call(rbind, changes) else data.frame()
  write.table(ch, OUTPUT_FILE_CURATION_CHANGES, sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(do.call(rbind, summary_rows), OUTPUT_FILE_CURATION_SUMMARY,
              sep = "\t", quote = FALSE, row.names = FALSE)
  mcra_log(sprintf("  dada2 training set: %d sequences → %s", length(ts$seq), OUTPUT_FILE_CURATED_TRAIN))
  mcra_log(sprintf("  dada2 species set : %d sequences → %s", length(ss$id), OUTPUT_FILE_CURATED_SPECIES))
  if (nrow(ch)) for (s in unique(ch$stage))
    mcra_log(sprintf("  label changes %-24s %d", s, sum(ch$stage == s)))

  list(lineages = setNames(lapply(seq_len(nrow(lab_final)), function(i) lab_final[i, TAX_LEVELS]),
                           db$safe_id))
}

# =============================================================================
# STEP 5d — Template sequence merging (fill taxa gaps in current DB)
# =============================================================================
#
# Loads ALL sequences from mcrAtemplate.fasta and BLASTs them against the
# current NCBI-derived candidate pool.  Any template sequence with no
# 100%-identical (full-length) match in the current pool is returned as a
# pre_trim_list-compatible entry, annotated with its tax4mcrA.taxonomy
# lineage.  These unique sequences are subsequently verified with taxize and
# passed through the same length / duplicate filters in Step 6, ensuring that
# no taxon present in the original Yang et al. template is silently lost.

find_unique_template_seqs <- function(template_fasta, template_taxonomy,
                                      current_candidates) {
  if (!file.exists(template_fasta)) {
    mcra_log(sprintf("  Template FASTA '%s' not found — merging skipped", template_fasta), "WARN")
    return(list())
  }

  tmpl_raw  <- read_fasta(template_fasta)
  tmpl_seqs <- setNames(
    toupper(gsub("[^A-Za-z]", "", tmpl_raw)),
    names(tmpl_raw)
  )
  mcra_log(sprintf("  Template sequences loaded: %d", length(tmpl_seqs)))

  # Safe IDs for template sequences (needed regardless of whether BLAST runs)
  safe_tmpl    <- gsub("[|;,\\[\\](){}<>\\\\/\\s]", "_", names(tmpl_seqs), perl = TRUE)
  safe_to_orig <- setNames(names(tmpl_seqs), safe_tmpl)

  # ── Determine unique template sequences ─────────────────────────────────────
  # If the current candidate pool is empty (e.g. HMMER rejected everything),
  # there is nothing to compare against — every template sequence is unique.
  # Skip BLAST entirely in that case.
  unique_orig <- if (length(current_candidates) == 0L) {

    mcra_log(paste0("  Current candidate pool is empty — no BLAST needed;",
                    " all template sequences treated as unique."), "WARN")
    names(tmpl_seqs)   # every template sequence is unique by definition

  } else {

    tmp_dir    <- tempdir()
    curr_fa    <- file.path(tmp_dir, "merge_curr.fa")
    tmpl_fa    <- file.path(tmp_dir, "merge_tmpl.fa")
    curr_db    <- file.path(tmp_dir, "merge_curr_db")
    result_tsv <- file.path(tmp_dir, "merge_tmpl_vs_curr.tsv")

    safe_curr <- gsub("[|;,\\[\\](){}<>\\\\/\\s]", "_", names(current_candidates), perl = TRUE)
    write_fasta(current_candidates, curr_fa, headers = safe_curr)
    write_fasta(tmpl_seqs,          tmpl_fa, headers = safe_tmpl)

    if (!build_blast_db(curr_fa, curr_db, label = "current candidates")) {
      mcra_log("  Falling back: treating all template sequences as unique.", "WARN")
      names(tmpl_seqs)
    } else {
      # blastn: template queries vs current DB at 100 % identity
      cmd_blast <- sprintf(
        paste("%s -query %s -db %s -out %s",
              "-outfmt '6 qseqid sseqid pident length qlen'",
              "-perc_identity 100 -max_target_seqs 1 -num_threads %d -dust no"),
        BLASTN_PATH, shQuote(tmpl_fa), shQuote(curr_db),
        shQuote(result_tsv), BLAST_THREADS
      )
      system(cmd_blast, ignore.stdout = TRUE, ignore.stderr = TRUE)

      matched_safe <- character(0)
      if (file.exists(result_tsv) && file.info(result_tsv)$size > 0L) {
        hits <- tryCatch(
          read.table(result_tsv, sep = "\t",
                     col.names = c("qseqid", "sseqid", "pident", "length", "qlen"),
                     stringsAsFactors = FALSE),
          error = function(e) data.frame()
        )
        if (nrow(hits) > 0L) {
          hits <- hits[hits$pident >= 100 & hits$length >= hits$qlen, ]
          matched_safe <- unique(hits$qseqid)
        }
      }
      safe_to_orig[setdiff(safe_tmpl, matched_safe)]   # orig accessions of unique seqs
    }
  }

  mcra_log(sprintf(
    "  Template: %d total | %d already in current (100%% identical) | %d unique → to be added",
    length(tmpl_seqs),
    length(tmpl_seqs) - length(unique_orig),
    length(unique_orig)
  ))

  if (length(unique_orig) == 0L) return(list())

  # ── Build pre_trim_list-compatible entries ──────────────────────────────────
  entries   <- list()
  n_no_tax  <- 0L
  n_too_short <- 0L

  for (orig_acc in unique_orig) {
    seq_str <- tmpl_seqs[[orig_acc]]
    if (is.null(seq_str) || is.na(seq_str) || nchar(seq_str) < MIN_LENGTH) {
      n_too_short <- n_too_short + 1L
      next
    }
    lineage <- template_taxonomy[[orig_acc]]
    if (is.null(lineage)) {
      n_no_tax <- n_no_tax + 1L
      next
    }
    entries[[orig_acc]] <- list(
      accession     = orig_acc,
      header        = format_dada2_header(lineage, accession = orig_acc),
      seq_str       = seq_str,
      lineage       = lineage,
      yang_flag     = FALSE,
      template_flag = TRUE    # distinguishes these as template-sourced additions
    )
  }

  mcra_log(sprintf(
    "  Unique template entries prepared: %d  (skipped: %d no taxonomy, %d too short)",
    length(entries), n_no_tax, n_too_short
  ))
  entries
}


# =============================================================================
# MAIN PIPELINE
# =============================================================================

main <- function() {
  mcra_log(paste(rep("=", 60), collapse = ""))
  mcra_log("mcrA DADA2 Database Builder (R version) — starting")
  mcra_log(paste(rep("=", 60), collapse = ""))

  # Set NCBI credentials
  options(entrez_email = ENTREZ_EMAIL)
  if (nchar(ENTREZ_API_KEY) > 0) {
    mcra_log("NCBI API key detected — 10 req/sec rate limit")
  } else {
    mcra_log("No NCBI API key — 3 req/sec rate limit")
  }

  # ── Step 0: Cluster seed template at 95 % ───────────────────────────────────
  mcra_log(sprintf("\n[Step 0] Clustering seed template at %.0f%% identity (cd-hit-est)",
                   CDHIT_SEED_IDENTITY * 100))
  cluster_seeds(SEED_TEMPLATE, SEED_FILE)

  # ── Step 1: Load seeds ───────────────────────────────────────────────────────
  mcra_log("\n[Step 1] Loading seed sequences")
  seeds <- load_seeds(SEED_FILE)

  # ── Step 2: BLAST seeds ──────────────────────────────────────────────────────
  mcra_log(sprintf("\n[Step 2] BLASTing seeds against NCBI '%s'", BLAST_DB))
  if (nchar(ENTREZ_QUERY) > 0) mcra_log(sprintf("  Entrez filter: %s", ENTREZ_QUERY))
  accession_list <- blast_all_seeds(seeds)
  if (length(accession_list) == 0) stop("No BLAST hits found.")

  # ── Step 3: Fetch sequences ──────────────────────────────────────────────────
  mcra_log("\n[Step 3] Fetching sequences from NCBI")

  if (file.exists(SEQ_CACHE_FILE)) {
    # Level-2 cache: load fully-parsed seq_df and skip fetch+parse entirely
    mcra_log(sprintf("  Found parsed seq_df cache '%s' — loading (skipping fetch+parse).",
                     SEQ_CACHE_FILE))
    mcra_log("  Delete this file to force a fresh fetch.")
    seq_df <- readRDS(SEQ_CACHE_FILE)
    mcra_log(sprintf("  Loaded %d rows from seq_df cache.", nrow(seq_df)))
  } else {
    # Level-1 cache: per-batch GenBank files in GENBANK_CACHE_DIR
    gb_batches <- fetch_sequences(accession_list)
    seq_df     <- parse_genbank_records(gb_batches)
    if (nrow(seq_df) == 0) stop("No sequences retrieved from NCBI.")
    mcra_log(sprintf("  Parsed %d GenBank records", nrow(seq_df)))

    # Persist parsed result so downstream failures don't require re-fetching
    tryCatch({
      saveRDS(seq_df, SEQ_CACHE_FILE)
      mcra_log(sprintf("  Saved parsed seq_df to cache '%s'", SEQ_CACHE_FILE))
    }, error = function(e) {
      mcra_log(sprintf("  Could not write seq_df cache: %s", e$message), "WARN")
    })
  }

  # ── Step 4: Fetch and correct taxonomy ──────────────────────────────────────
  mcra_log("\n[Step 4] Fetching taxonomy from NCBI")

  if (file.exists(TAX_CACHE_FILE)) {
    mcra_log(sprintf("  Found taxonomy cache '%s' — loading (skipping fetch).", TAX_CACHE_FILE))
    mcra_log("  Delete this file to force a fresh taxonomy fetch.")
    taxid_lineage <- readRDS(TAX_CACHE_FILE)
    mcra_log(sprintf("  Loaded %d taxon entries from cache.", length(taxid_lineage)))
  } else {
    taxid_lineage <- fetch_taxonomy_batch(seq_df$taxid)

    if (CORRECT_TAXONOMY) {
      mcra_log("  Verifying taxonomy names against current NCBI (taxize)...")
      n_corrected <- 0
      for (tid in names(taxid_lineage)) {
        old_sp  <- taxid_lineage[[tid]]["species"]
        new_lin <- correct_lineage_names(taxid_lineage[[tid]])
        taxid_lineage[[tid]] <- new_lin
        if (new_lin["species"] != old_sp) n_corrected <- n_corrected + 1
      }
      mcra_log(sprintf("  Taxonomy correction: %d name(s) updated", n_corrected))
    }

    tryCatch({
      saveRDS(taxid_lineage, TAX_CACHE_FILE)
      mcra_log(sprintf("  Saved taxonomy to cache '%s'", TAX_CACHE_FILE))
    }, error = function(e) {
      mcra_log(sprintf("  Could not write taxonomy cache: %s", e$message), "WARN")
    })
  }

  # ── Step 5: HMMER gene-coverage filter ──────────────────────────────────────
  hmm_trimmed  <- list()
  n_hmm_reject <- 0

  if (!is.null(HMM_PROFILE) && file.exists(HMM_PROFILE)) {
    if (file.exists(HMMER_CACHE_FILE)) {
      mcra_log(sprintf(
        "\n[Step 5] HMMER — loading checkpoint '%s' (skipping hmmsearch).",
        HMMER_CACHE_FILE))
      mcra_log("  Delete this file to re-run HMMER with a different coverage threshold.")
      hmm_trimmed  <- readRDS(HMMER_CACHE_FILE)
      n_hmm_reject <- nrow(seq_df) - length(hmm_trimmed)
      mcra_log(sprintf("  Loaded %d HMM-trimmed sequences from checkpoint.",
                       length(hmm_trimmed)))
    } else {
      mcra_log(sprintf(
        "\n[Step 5] HMMER gene-coverage filter + trimming (min_coverage=%.0f%%)",
        HMM_MIN_COVERAGE * 100))
      hmm_result   <- run_hmmsearch_filter(seq_df, HMM_PROFILE)
      hmm_trimmed  <- if (!is.null(hmm_result)) hmm_result else list()
      n_hmm_reject <- nrow(seq_df) - length(hmm_trimmed)
      tryCatch({
        saveRDS(hmm_trimmed, HMMER_CACHE_FILE)
        mcra_log(sprintf("  Checkpoint saved → '%s'", HMMER_CACHE_FILE))
      }, error = function(e) {
        mcra_log(sprintf("  Could not write HMMER checkpoint: %s", e$message), "WARN")
      })
    }
  } else {
    mcra_log("\n[Step 5] HMMER skipped (HMM_PROFILE is NULL or file not found)")
    if (EXTRACT_COMPLETE_GENE) {
      mcra_log("  Falling back to GenBank CDS extraction")
    } else {
      mcra_log("  Using raw fetched sequences")
    }
  }

  # ── Step 5b/5c/5d replaced — see Step 6 merge below ─────────────────────────

  # ── Step 6: Merge, filter, deduplicate, build DB entries ────────────────────
  # Sources:
  #   1. BLAST-retrieved sequences (hmm_trimmed from NCBI)
  #   2. yang_2014_taxize.fasta  — Yang 2014 with NCBI-enriched lineages
  #   3. mcrA_ncbi_genome_db.fasta — NCBI genome annotation-based sequences
  #
  # Deduplication: for exact-duplicate sequences keep the entry with the most
  # comprehensive lineage (fewest Unclassified / uncultured rank values).
  # ─────────────────────────────────────────────────────────────────────────────
  mcra_log("\n[Step 6] Merging databases, filtering, deduplicating, building DB entries")

  # ── Helper: parse DADA2-format FASTA into pre_trim_list entries ──────────────
  parse_dada2_fasta_to_pool <- function(fasta_path, source_label) {
    if (!file.exists(fasta_path)) {
      mcra_log(sprintf("  [%s] not found — skipped", fasta_path), "WARN"); return(list())
    }
    lines   <- readLines(fasta_path)
    hdr_idx <- grep("^>", lines)
    pool    <- list()
    for (i in seq_along(hdr_idx)) {
      hdr <- sub("^>\\s*", "", lines[hdr_idx[i]])
      seq_lines <- lines[(hdr_idx[i]+1L):(if (i < length(hdr_idx)) hdr_idx[i+1L]-1L else length(lines))]
      seq <- toupper(gsub("[^ACGTNacgtn]", "", paste(seq_lines, collapse="")))
      if (!nzchar(seq)) next
      parts <- strsplit(hdr, ";", fixed=TRUE)[[1L]]
      n     <- length(parts)
      acc   <- if (n >= 8L) trimws(parts[8L]) else paste0(source_label, "_", i)
      lin   <- setNames(rep("Unclassified", length(TAX_LEVELS)), TAX_LEVELS)
      rank_map <- c(Domain="kingdom",Phylum="phylum",Class="class",
                    Order="order",Family="family",Genus="genus",Species="species")
      for (j in seq_len(min(n-1L,7L))) {
        val <- trimws(parts[j])
        if (nzchar(val) && val != "Unclassified") lin[rank_map[names(rank_map)[j]]] <- val
      }
      header <- format_dada2_header(lin, accession=acc)
      pool[[acc]] <- list(accession=acc, header=header, seq_str=seq, lineage=lin,
                          yang_flag=FALSE, template_flag=FALSE, source=source_label)
    }
    mcra_log(sprintf("  [%s] loaded: %d sequences", source_label, length(pool)))
    pool
  }

  lineage_score <- function(lin) {
    if (is.null(lin)) return(0L)
    UNRES <- "^(Unclassified|unclassified|uncultured|NA|)$|^uncultured|_unclassified$"
    sum(!is.na(lin) & nzchar(lin) & !grepl(UNRES, lin))
  }

  blast_pool <- list()
  for (i in seq_len(nrow(seq_df))) {
    acc <- seq_df$accession[i]
    seq_str <- if (!is.null(HMM_PROFILE) && file.exists(HMM_PROFILE)) hmm_trimmed[[acc]]
               else if (EXTRACT_COMPLETE_GENE) seq_df$cds_seq[i] else seq_df$sequence[i]
    if (is.null(seq_str)||is.na(seq_str)||nchar(seq_str)<MIN_LENGTH) next
    taxid   <- seq_df$taxid[i]
    lineage <- if (!is.na(taxid)) taxid_lineage[[taxid]] else NULL
    if (is.null(lineage)) lineage <- setNames(rep("Unclassified",length(TAX_LEVELS)),TAX_LEVELS)
    header  <- format_dada2_header(lineage, accession=acc)
    blast_pool[[acc]] <- list(accession=acc,header=header,seq_str=seq_str,lineage=lineage,
                               yang_flag=FALSE,template_flag=FALSE,source="BLAST")
  }
  mcra_log(sprintf("  BLAST pool: %d sequences", length(blast_pool)))

  yang_pool       <- parse_dada2_fasta_to_pool(YANG_TAXIZE_FILE,  "Yang2014")
  ncbi_annot_pool <- parse_dada2_fasta_to_pool(NCBI_ANNOT_FILE,   "NCBIannot")
  all_candidates  <- c(blast_pool, yang_pool, ncbi_annot_pool)
  mcra_log(sprintf("  Combined pool before dedup: %d sequences", length(all_candidates)))

  seq_to_accs <- list()
  for (acc in names(all_candidates)) {
    s <- toupper(gsub("[^ACGTN]","",all_candidates[[acc]]$seq_str))
    seq_to_accs[[s]] <- c(seq_to_accs[[s]], acc)
  }

  pre_trim_list <- list(); acc_corr <- list()
  n_kept <- 0L; n_dropped <- 0L; n_lfilt <- 0L; n_notax <- 0L

  for (s in names(seq_to_accs)) {
    group  <- seq_to_accs[[s]]
    scores <- vapply(group, function(a) lineage_score(all_candidates[[a]]$lineage), integer(1L))
    best   <- group[which.max(scores)]
    entry  <- all_candidates[[best]]
    filt   <- passes_filters(entry$seq_str)
    if (!filt$ok) {
      n_lfilt <- n_lfilt+1L
      acc_corr[[length(acc_corr)+1L]] <- list(accession=entry$accession,dada2_header="",
        seq_length=nchar(entry$seq_str),template_annotated=FALSE,yang_annotated=entry$yang_flag,
        in_final_db=FALSE,status="length_filtered"); next
    }
    if (is.null(entry$lineage)) {
      n_notax <- n_notax+1L
      acc_corr[[length(acc_corr)+1L]] <- list(accession=entry$accession,dada2_header="",
        seq_length=nchar(entry$seq_str),template_annotated=FALSE,yang_annotated=entry$yang_flag,
        in_final_db=FALSE,status="no_taxonomy"); next
    }
    pre_trim_list[[entry$accession]] <- entry; n_kept <- n_kept+1L
    for (dup in setdiff(group,best)) {
      n_dropped <- n_dropped+1L
      acc_corr[[length(acc_corr)+1L]] <- list(accession=all_candidates[[dup]]$accession,
        dada2_header="",seq_length=nchar(entry$seq_str),template_annotated=FALSE,
        yang_annotated=all_candidates[[dup]]$yang_flag,in_final_db=FALSE,
        status="duplicate_lower_score")
    }
  }
  mcra_log(sprintf("  After dedup: %d kept, %d dropped, %d length-filtered, %d no-taxonomy",
                   n_kept,n_dropped,n_lfilt,n_notax))
  mcra_log(sprintf("  Total pre-trimming entries: %d", length(pre_trim_list)))

  # ── Step 7: Alignment-based gap trimming ────────────────────────────────────
  mcra_log("\n[Step 7] Alignment-based gap trimming (MAFFT)")

  seq_pool <- setNames(
    sapply(pre_trim_list, function(e) e$seq_str),
    names(pre_trim_list)
  )
  trim_result <- trim_by_gap_profile(seq_pool)
  trimmed_seqs <- trim_result$kept
  gap_removed  <- trim_result$removed

  # ── Step 8: Write final DADA2 database ──────────────────────────────────────
  mcra_log(sprintf("\n[Step 8] Writing final DADA2 database → %s", OUTPUT_FILE))

  fh      <- file(OUTPUT_FILE, "w")
  written <- 0

  for (acc in names(trimmed_seqs)) {
    entry <- pre_trim_list[[acc]]
    if (is.null(entry)) next

    trimmed_seq <- trimmed_seqs[[acc]]
    filt <- passes_filters(trimmed_seq)
    if (!filt$ok) {
      gap_removed[[acc]] <- sprintf("post-trim filter: %s", filt$reason)
      next
    }

    cat(sprintf("%s\n%s\n", entry$header, trimmed_seq), file = fh)
    written <- written + 1
    acc_corr[[length(acc_corr)+1]] <- list(
      accession          = acc,
      dada2_header       = sub("^>", "", entry$header),
      seq_length         = nchar(trimmed_seq),
      template_annotated = isTRUE(entry$template_flag),
      yang_annotated     = isTRUE(entry$yang_flag),
      in_final_db        = TRUE,
      status             = "accepted"
    )
  }
  close(fh)

  # Record gap-removed sequences
  for (acc in names(gap_removed)) {
    entry <- pre_trim_list[[acc]]
    acc_corr[[length(acc_corr)+1]] <- list(
      accession          = acc,
      dada2_header       = "",
      seq_length         = if (!is.null(entry)) nchar(entry$seq_str) else 0,
      template_annotated = if (!is.null(entry)) isTRUE(entry$template_flag) else FALSE,
      yang_annotated     = if (!is.null(entry)) isTRUE(entry$yang_flag) else FALSE,
      in_final_db        = FALSE,
      status             = "gap_trimmed"
    )
  }

  mcra_log(sprintf("  Final database: %d sequences → %s", written, OUTPUT_FILE))

  # ── Step 9: Phylogenetic taxonomy curation ──────────────────────────────────
  # Build a tree from the accepted sequences; label internal nodes where all
  # known descendants agree on a taxon AND the max patristic distance is within
  # the rank threshold; propagate those labels to "Unclassified" tips;
  # re-verify with taxize; write curated FASTA + Newick tree for external review.

  curated_pre_trim <- pre_trim_list   # will be updated in-place if curation runs

  if (PHYLO_CURATE_TAXONOMY && written > 0L) {
    mcra_log("\n[Step 9] Phylogenetic taxonomy curation")

    # ── Collect final sequences and lineages (accepted by Step 8) ────────────
    final_accs <- intersect(names(trimmed_seqs), names(pre_trim_list))
    final_seqs <- setNames(
      sapply(final_accs, function(a) trimmed_seqs[[a]]),
      final_accs
    )
    final_lins <- setNames(
      lapply(final_accs, function(a) pre_trim_list[[a]]$lineage),
      final_accs
    )

    # Safe IDs for tree tip labels (blastn / MAFFT header constraints)
    safe_ids   <- gsub("[|;,\\[\\](){}<>\\\\/\\s]", "_", final_accs, perl = TRUE)
    seen_safe2 <- character(0)
    for (i in seq_along(safe_ids)) {
      base <- safe_ids[i]; k <- 0L
      while (safe_ids[i] %in% seen_safe2) {
        k <- k + 1L; safe_ids[i] <- sprintf("%s_%d", base, k)
      }
      seen_safe2 <- c(seen_safe2, safe_ids[i])
    }
    safe_to_orig2 <- setNames(final_accs, safe_ids)
    safe_pool2    <- setNames(unname(final_seqs), safe_ids)
    safe_lins2    <- setNames(final_lins,         safe_ids)

    # ── Build tree ────────────────────────────────────────────────────────────
    tree_ok <- build_mcra_tree(safe_pool2, TREE_FILE)

    if (tree_ok) {
      phy_tree <- tryCatch(ape::read.tree(TREE_FILE), error = function(e) NULL)

      if (!is.null(phy_tree)) {
        mcra_log(sprintf("  Tree: %d tips, %d internal nodes",
                         length(phy_tree$tip.label), phy_tree$Nnode))

        # ── 9a cleaning → 9b tree curation → 9c fill ─────────────────────────
        curate_res <- tryCatch(
          refine_and_curate_step9(final_accs, safe_pool2, safe_lins2, phy_tree),
          error = function(e) { mcra_log(sprintf("  Curation failed: %s", e$message), "WARN"); NULL }
        )

        if (!is.null(curate_res)) {
          # ── Write curated FASTA ───────────────────────────────────────────
          mcra_log(sprintf("  Writing curated FASTA → %s", OUTPUT_FILE_CURATED))
          fh_cur      <- file(OUTPUT_FILE_CURATED, "w")
          written_cur <- 0L
          for (safe_id in names(safe_pool2)) {
            orig_acc <- safe_to_orig2[safe_id]
            entry    <- pre_trim_list[[orig_acc]]
            seq_str  <- trimmed_seqs[[orig_acc]]
            if (is.null(entry) || is.null(seq_str)) next

            lin <- curate_res$lineages[[safe_id]]
            if (is.null(lin)) lin <- entry$lineage
            hdr <- format_dada2_header(lin, accession = orig_acc)
            cat(sprintf("%s\n%s\n", hdr, seq_str), file = fh_cur)

            # Update curated_pre_trim so Step 10 correspondence reflects curation
            curated_pre_trim[[orig_acc]]$lineage <- lin
            curated_pre_trim[[orig_acc]]$header  <- hdr

            written_cur <- written_cur + 1L
          }
          close(fh_cur)
          mcra_log(sprintf("  Curated database: %d sequences → %s",
                           written_cur, OUTPUT_FILE_CURATED))
          mcra_log(sprintf("  Tree saved for review → %s  (open in FigTree / iTOL)",
                           TREE_FILE))
        } else {
          mcra_log("  Curation failed — curated FASTA not written", "WARN")
        }
      } else {
        mcra_log("  Could not read tree file — curation skipped", "WARN")
      }
    } else {
      mcra_log("  Tree build failed — phylo curation skipped", "WARN")
    }
  } else if (!PHYLO_CURATE_TAXONOMY) {
    mcra_log("\n[Step 9] Phylogenetic taxonomy curation skipped (PHYLO_CURATE_TAXONOMY = FALSE)")
  }

  # ── Step 10: Write correspondence file ──────────────────────────────────────
  mcra_log(sprintf("\n[Step 10] Writing correspondence file → %s", OUTPUT_CORRESPONDENCE))
  write_correspondence_tsv(acc_corr, OUTPUT_CORRESPONDENCE)

  # ── Summary ──────────────────────────────────────────────────────────────────
  n_tmpl <- sum(sapply(acc_corr, function(r) isTRUE(r$template_annotated) && isTRUE(r$in_final_db)))
  n_yang <- sum(sapply(acc_corr, function(r) isTRUE(r$yang_annotated) && isTRUE(r$in_final_db)))
  mcra_log(paste(rep("=", 60), collapse = ""))
  mcra_log("DONE")
  mcra_log(sprintf("  Seeds clustered (95 %%)  : %s", SEED_FILE))
  mcra_log(sprintf("  Final database           : %s  (%d sequences)", OUTPUT_FILE, written))
  if (PHYLO_CURATE_TAXONOMY)
    mcra_log(sprintf("  Curated database         : %s", OUTPUT_FILE_CURATED))
  mcra_log(sprintf("  Tree (for review)        : %s", TREE_FILE))
  mcra_log(sprintf("  Correspondence file      : %s", OUTPUT_CORRESPONDENCE))
  mcra_log(sprintf("  Template-annotated seqs  : %d", n_tmpl))
  mcra_log(sprintf("  Yang-annotated seqs      : %d", n_yang))
  mcra_log(sprintf("  NCBI-annotated seqs      : %d", written - n_tmpl - n_yang))
  mcra_log(sprintf("  Removed by gap trimming  : %d", length(gap_removed)))
  mcra_log(sprintf("  Run log                  : %s", LOG_FILE))
  mcra_log(paste(rep("=", 60), collapse = ""))
  mcra_log("\nTo use in DADA2 (pre-curation):")
  mcra_log(sprintf('  taxa <- assignTaxonomy(seqs, "%s", multithread = TRUE,', OUTPUT_FILE))
  mcra_log('                         taxLevels = c("Kingdom","Phylum","Class",')
  mcra_log('                                       "Order","Family","Genus","Species","Accession"))')
  if (PHYLO_CURATE_TAXONOMY) {
    mcra_log("\nTo use curated database in DADA2:")
    mcra_log(sprintf('  taxa <- assignTaxonomy(seqs, "%s", multithread = TRUE)', OUTPUT_FILE_CURATED_TRAIN))
    mcra_log(sprintf('  taxa <- addSpecies(taxa, "%s")', OUTPUT_FILE_CURATED_SPECIES))
  }
}

# Run
main()
