# update_mcrA_db — Updated *mcrA* Reference Databases for DADA2

Curated *mcrA* (methyl-coenzyme M reductase subunit alpha) reference databases
in dada2, Mothur, and QIIME2 format, targeted at the amplicon window defined by the primer set of
[Angel et al. (2012)](https://doi.org/10.1038/ismej.2011.141), together with
the R scripts used to build and evaluate them.

The databases supersede the Yang et al. (2014) *mcrA* reference set with
broader sequence coverage, HMM-verified gene boundaries, and updated NCBI
taxonomy.

---

## Repository layout

```
update_mcrA_db/
├── databases/
│   ├── mcrA_ncbi_nt_db/               # updated DB — amplicon-length mcrA (27,942 seqs)
│   │   ├── dada2/   mcrA_ncbi_nt_db.fasta
│   │   ├── qiime2/  mcrA_ncbi_nt_db_seqs.fasta
│   │   │            mcrA_ncbi_nt_db_taxonomy.tsv
│   │   └── mothur/  mcrA_ncbi_nt_db.fasta
│   │                mcrA_ncbi_nt_db.taxonomy
│   ├── mcrA_ncbi_nt_cur_db/           # updated DB with curated taxonomy  (27,942 seqs)
│   │   ├── dada2/   mcrA_ncbi_nt_cur_db.fasta
│   │   ├── qiime2/  mcrA_ncbi_nt_cur_db_seqs.fasta
│   │   │            mcrA_ncbi_nt_cur_db_taxonomy.tsv
│   │   └── mothur/  mcrA_ncbi_nt_cur_db.fasta
│   │                mcrA_ncbi_nt_cur_db.taxonomy
│   └── mcrA_ncbi_genome_db/           # NCBI genome annotation-based DB   (1,572 seqs)
│       ├── dada2/   mcrA_ncbi_genome_db.fasta
│       ├── qiime2/  mcrA_ncbi_genome_db_seqs.fasta
│       │            mcrA_ncbi_genome_db_taxonomy.tsv
│       └── mothur/  mcrA_ncbi_genome_db.fasta
│                    mcrA_ncbi_genome_db.taxonomy
├── scripts/
│   ├── build_mcrA_db_ncbi_annot.R     # Step 1 — NCBI genome DB + Yang 2014 DADA2 format
│   ├── calibrate_mcra_thresholds.R    # Step 2 — phylogenetic rank-threshold calibration
│   ├── build_mcrA_db.R                # Step 3 — amplicon DB (seed-and-expand from NCBI nt)
│   ├── compare_mcrA_databases2.R      # Step 4 — cross-database comparison
│   ├── convert_db_formats.R           # Step 5 — convert DADA2 databases to QIIME2 / Mothur
│   ├── run_virtualPCR.sh              # Step 6a — in silico PCR (single mismatch level)
│   ├── run_mismatch_sweep.sh          # Step 6b — sweep over 0–3 3'-end mismatches
│   ├── parse_virtualPCR.R             # Step 6c — parse single-run results → coverage tables
│   └── parse_mismatch_sweep.R         # Step 6d — parse sweep results → tables & plots
├── results/                           # Output of Step 6 scripts (created at run time)
│   ├── primer_coverage/               # Output of run_virtualPCR.sh
│   ├── primer_coverage_sweep_errors0/ # Output of mismatch sweep at 0 mismatches
│   ├── primer_coverage_sweep_errors1/ #                              1
│   ├── primer_coverage_sweep_errors2/ #                              2
│   ├── primer_coverage_sweep_errors3/ #                              3
│   ├── mismatch_sweep_summary.tsv     # Long-format cross-mismatch coverage table
│   ├── mismatch_sweep_wide.tsv        # Wide-format pivot
│   └── taxonomy/                      # Per-taxon TSVs and multi-panel PDF plot
├── software_versions.tsv              # R package and external tool versions
└── README.md
```

---

## Databases

| Database folder       | Label (manuscript) | Sequences | Description                                                                                                                                                                                        |
| --------------------- | :----------------: | --------: | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `mcrA_ncbi_nt_db`     | **mlas-MCRA**      |    27,942 | Amplicon-length mcrA sequences retrieved from NCBI nt via BLAST against the Yang et al. 2014 seed template; HMMER-verified (TIGR03256); taxonomy resolved with taxize against current NCBI records |
| `mcrA_ncbi_nt_cur_db` | **mlas-MCRA_CUR**  |    27,942 | Same sequences as `mcrA_ncbi_nt_db`; taxonomy additionally curated by phylogenetic inference — a GTR tree is built with FastTree; internal nodes are labelled at rank R with taxon X only when (i) all classified descendants agree unanimously on X, (ii) max intra-subtree patristic distance ≤ calibrated rank threshold (Step 2), and (iii) all tips carrying X in the full tree are contained within the subtree (strict monophyly); labels are propagated to "Unclassified" tips from the deepest qualifying ancestor; updated names are re-verified with taxize |
| `mcrA_ncbi_genome_db` | **MCRA_FULL**      |     1,572 | Full CDS-length mcrA sequences extracted directly from annotated NCBI genomes of known methanogenic archaea; all sequences verified against TIGR03256                                              |

Each database folder contains three format subfolders: `dada2/`, `qiime2/`, and `mothur/`
(see [Format conversion](#5-convert_db_formatsr) below).

> **Note — amplicon vs full-length sequences.** `mcrA_ncbi_nt_db` and
> `mcrA_ncbi_nt_cur_db` contain sequences in the amplicon window of the
> Angel et al. (2012) primers (~450–470 bp). `mcrA_ncbi_genome_db`
> contains full-length mcrA CDS sequences (~1,650 bp) and is most useful as a
> phylogenetic reference or for primer evaluation against complete genes.


---

## Scripts

Scripts are intended to be run **in the order listed below**. Each script
produces outputs that are consumed by subsequent steps.

### 1. `build_mcrA_db_ncbi_annot.R`

Builds the NCBI genome annotation-based mcrA database
(`mcrA_ncbi_genome_db.fasta`) and converts the Yang et al. (2014) mcrA
reference set into DADA2 format (`yang_2014_taxize.fasta`).

**Workflow:**

1. Queries the NCBI Datasets v2 API for all methanogenic archaeal genome
   accessions across the major methanogenic lineages (Methanobacteria,
   Methanococci, Methanomicrobia, Methanopyri, and ANME clades)
2. Fetches complete taxonomic lineages via NCBI Taxonomy (rentrez)
3. Downloads annotated CDS nucleotide FASTA per genome
4. Extracts mcrA-annotated CDS by header pattern matching
   (`[gene=mcrA]` or *methyl-coenzyme M reductase alpha*)
5. Translates each CDS and runs hmmsearch against `TIGR03256.hmm`
6. Retains sequences covering ≥ 90 % of the HMM profile length
7. Trims nucleotide sequences to HMM envelope coordinates
8. Deduplicates by exact trimmed nucleotide sequence
9. Formats DADA2 taxonomy headers and writes `mcrA_ncbi_genome_db.fasta`
10. Reformats the Yang et al. (2014) template with NCBI-resolved lineages and
    writes `yang_2014_taxize.fasta`

**Required input files:**

- `TIGR03256.hmm` — TIGRFAM HMM profile for methyl-coenzyme M reductase alpha
  (download from [TIGRFAM](https://tigrfams.jcvi.org/cgi-bin/HmmReportPage.cgi?acc=TIGR03256)
  or via HMMER's hmmfetch from the TIGRFAM database)
- `mcrAtemplate.fasta` — Yang et al. (2014) mcrA reference sequences (FASTA)
- `tax4mcrA.taxonomy` — Yang et al. (2014) taxonomy file

**R packages:** `httr`, `jsonlite`, `rentrez`, `Biostrings`, `dplyr`, `stringr`

**External tools:** `hmmer` (hmmsearch, hmmpress), `ncbi-datasets-cli` (optional)

---

### 2. `calibrate_mcra_thresholds.R`

Empirically calibrates phylogenetic patristic distance thresholds for each
taxonomic rank (species → phylum) using the Youden-J ROC method of [Valles-Colomer et al. (2023)](https://doi.org/10.1038/s41586-022-05620-1), applied here to rank-level patristic distances following the framework of [Yarza et al. (2014)](https://doi.org/10.1038/nrmicro3330).

**Workflow:**

1. Loads sequences and taxonomy of the Yang et al. (2014) database: `mcrAtemplate.fasta` + `tax4mcrA.taxonomy`
2. Aligns sequences with MAFFT; builds a GTR tree with FastTree (NJ fallback)
3. Computes the full patristic distance matrix
4. For every sequence pair, records the deepest rank at which the two
   sequences differ taxonomically
5. For each rank, performs a binary classification (inter-rank vs intra-rank
   pairs) and finds the distance cutoff *d\** that maximises the Youden index:
   *J(d)* = sensitivity(*d*) + specificity(*d*) − 1
6. Reports *d\** per rank with 95 % bootstrap confidence intervals (n = 999)

The resulting thresholds are hardcoded as constants into `build_mcrA_db.R`
(Step 3) and used in Step 9 to constrain phylogenetic taxonomy curation: an
internal tree node receives a rank label only when the maximum patristic distance
across all its known-taxonomy descendants does not exceed the rank threshold.
After re-running this script, update the `PHYLO_RANK_THRESHOLDS` constants in
`build_mcrA_db.R` manually.

**Required input files:**

- `mcrAtemplate.fasta` — Yang et al. (2014) mcrA reference sequences
- `tax4mcrA.taxonomy` — Yang et al. (2014) taxonomy file

**R packages:** `ape`, `phangorn`, `ggplot2`, `dplyr`, `tidyr`, `pROC`,
`gridExtra`

**External tools:** `mafft`, `FastTree`

**Outputs:** `mcra_rank_thresholds.tsv`, `mcra_rank_thresholds.pdf`

---

### 3. `build_mcrA_db.R`

Builds two versions of the amplicon-length mcrA database using a
seed-and-expand strategy: BLAST seeds against NCBI nt, fetch full sequences
and GenBank records, apply HMM coverage filters, transfer curated taxonomy,
and trim to the primer-amplified region.

**Workflow:**

0. Clusters the seed template at 95 % identity (CD-HIT-EST)
1. BLASTs web-seed sequences against NCBI nt (NCBI BLAST API via rentrez / httr)
2. Fetches full sequences and GenBank records (rentrez)
3. Retrieves and corrects taxonomy with taxize against current NCBI records
4. Applies HMMER gene-coverage filter and trims to HMM envelope coordinates
5. Merges BLAST-retrieved sequences with `yang_2014_taxize.fasta` and
   `mcrA_ncbi_genome_db.fasta`; deduplicates exact-sequence duplicates by
   retaining the entry with the most complete lineage; appends NCBI accession
   to the Species field
6. Performs MAFFT-based gap trimming of the final alignment
7. Writes `mcrA_ncbi_nt_db.fasta` (NCBI taxonomy)
8. Writes accession correspondence TSV
9. **Phylogenetic taxonomy curation** (→ `mcrA_ncbi_nt_cur_db.fasta`): builds a
   GTR tree with FastTree; computes the full patristic distance matrix
   (`ape::cophenetic.phylo`); labels each internal node at rank R with taxon X
   only if all three conditions hold simultaneously: (i) all classified
   (non-"Unclassified") tips within the subtree agree unanimously on taxon X at
   rank R; (ii) the maximum intra-subtree patristic distance is ≤ the calibrated
   rank threshold (genus/family/order/class — hardcoded from
   `calibrate_mcra_thresholds.R` output); and (iii) every tip in the **full
   tree** carrying taxon X at rank R is contained within that subtree (strict
   monophyly — paraphyletic nodes are rejected); for each "Unclassified" tip,
   walks tip-to-root and assigns ranks from the deepest ancestor that carries a
   label; coarser ranks still unresolved after that assignment are filled by
   unanimous consensus of all known sequences in the anchor's subtree; re-
   verifies newly assigned names with taxize; fills any remaining Unclassified
   higher ranks via taxize lineage lookup

**Required input files (produced by Steps 1 and 2):**

- `yang_2014_taxize.fasta` — Yang et al. 2014 in DADA2 format (Step 1)
- `mcrA_ncbi_genome_db.fasta` — NCBI genome annotation-based database (Step 1)
- `TIGR03256.hmm` — TIGRFAM HMM profile (same as Step 1)

> **Note:** the patristic-distance rank thresholds from Step 2
> (`mcra_rank_thresholds.tsv`) are **not read at runtime**; they are hardcoded
> as constants (`PHYLO_RANK_THRESHOLDS`) near the top of `build_mcrA_db.R`.
> After re-running Step 2, update those constants manually.

**R packages:** `rentrez`, `taxize`, `Biostrings`, `stringr`, `dplyr`,
`readr`, `httr`, `xml2`, `ape`, `phangorn`

**External tools:** `mafft`, `FastTree`, `hmmer` (hmmsearch), `blast` (blastn),
`cd-hit-est`

---

### 4. `compare_mcrA_databases2.R`

Compares the four mcrA databases across multiple metrics to evaluate coverage,
taxonomic resolution, and sequence composition.

**Analyses performed:**

- Sequence counts and lengths per database
- Taxonomic richness at each rank (Domain → Species)
- Rank resolution rates (% sequences with resolved taxonomy per rank)
- Venn diagram of shared sequences between databases
- Bray–Curtis dissimilarity / UniFrac-style overlap at genus and family level
- Pairwise heatmaps of shared genera/families
- Stacked bar charts of relative composition per rank

**Required input files (produced by Steps 1–3):**

- `yang_2014_taxize.fasta`
- `mcrA_ncbi_nt_db.fasta`
- `mcrA_ncbi_nt_cur_db.fasta`
- `mcrA_ncbi_genome_db.fasta`

**R packages:** `tidyverse`, `ggplot2`, `patchwork`, `scales`, `ggvenn`,
`pheatmap`, `vegan`

> **Note:** this script uses `rstudioapi::getSourceEditorContext()$path` to
> resolve its working directory and must be run interactively from within
> RStudio (via **Source**). It is not compatible with `Rscript` on the command
> line.

---

### 5. `convert_db_formats.R`

Converts each DADA2-format database into QIIME2 and Mothur formats and
organises the outputs in per-database format subfolders.

**Workflow:**

1. Scans the `databases/` directory for `.fasta` files
2. Parses every DADA2 header (`>Domain;…;Species;Accession`) and extracts the
   accession ID and seven-level taxonomy string
3. Disambiguates duplicate accession IDs (e.g. multiple CDS from the same
   genome assembly) by appending `_1`, `_2`, … suffixes
4. Writes three output subfolders per database:
   - **`dada2/`** — unmodified copy of the original DADA2 FASTA
   - **`qiime2/`** — accession-header FASTA + taxonomy TSV
     (`Feature ID ⇥ k__…;p__…;c__…;o__…;f__…;g__…;s__…`)
   - **`mothur/`** — accession-header FASTA + `.taxonomy` file
     (`SeqID ⇥ Domain;…;Species;` with trailing semicolon and spaces replaced
     by underscores)

**Usage:**

```r
# From within R (set working directory to the repo root first):
source("scripts/convert_db_formats.R")

# Or from the terminal:
Rscript scripts/convert_db_formats.R databases/
```

**Required input:** DADA2-format `.fasta` files in the target directory
(produced by Steps 1 and 3).

**No R packages required** — uses base R only.

---

### 6. In silico PCR — primer coverage assessment

Evaluates the coverage of the **mlas-mod-F / mcrA-rev-R** primer pair against
the three databases using [virtualPCR](https://doi.org/10.3389/fbinf.2024.1464197)
in silico PCR, then generates per-database and per-taxon coverage tables and plots.

> **Prerequisite — unzip databases before running.**
> The FASTA files are distributed as zip archives.  Unzip each one before running
> Step 6 scripts:
> ```bash
> cd databases
> for z in *.zip; do unzip "$z"; done
> cd ..
> ```

**Primers**

| Name | Sequence (5'→3') | Note |
|---|---|---|
| mlas-mod-F | `GGYGGTGTMGGDTTCACMCARTA` | Forward |
| mcrA-rev-R | `BGCGTAGTTVGGRTAGT` | Reverse; first 7 bases (CGTTCAT) of the published primer excluded |

#### Step 6a — Single in silico PCR run (`run_virtualPCR.sh`)

Runs virtualPCR with **1 allowed 3'-end mismatch** against all three databases.
Results are written to `results/primer_coverage/`.

```bash
bash scripts/run_virtualPCR.sh
```

Key parameters (edit at the top of the script):

| Variable | Default | Meaning |
|---|---|---|
| `N3_ERRORS` | `1` | Allowed mismatches at the 3' end of each primer |
| `MIN_LEN` | `100` | Minimum expected amplicon length (bp) |
| `MAX_LEN` | `900` | Maximum expected amplicon length (bp) |

Output per database: `results/primer_coverage/<db_name>/config.txt` and
`results/primer_coverage/<db_name>/report.out`.

#### Step 6b — Mismatch sweep (`run_mismatch_sweep.sh`)

Runs virtualPCR four times, varying `number3errors` from 0 to 3.
Results land in `results/primer_coverage_sweep_errors0/` through
`results/primer_coverage_sweep_errors3/`.

```bash
bash scripts/run_mismatch_sweep.sh            # runs all four levels (0–3)
bash scripts/run_mismatch_sweep.sh 0 2        # run only specific levels
```

#### Step 6c — Parse single-run results (`parse_virtualPCR.R`)

Parses `results/primer_coverage/` and writes:

- `results/primer_coverage/coverage_summary.tsv` — per-database overall coverage
- `results/primer_coverage/taxonomy/taxonomy_coverage_<db>.tsv` — per-taxon
  coverage (domain → genus) for each database
- `results/primer_coverage/taxonomy/plots/coverage_<level>.pdf` — line plots

```bash
Rscript scripts/parse_virtualPCR.R
```

#### Step 6d — Parse sweep results and generate plots (`parse_mismatch_sweep.R`)

Parses all `results/primer_coverage_sweep_errorsN/` directories and writes:

- `results/mismatch_sweep_summary.tsv` — long-format table (database × mismatches × coverage)
- `results/mismatch_sweep_wide.tsv` — wide-format pivot
- `results/taxonomy/sweep_taxonomy_<level>.tsv` — per-taxon coverage per mismatch level
- `results/taxonomy/sweep_taxonomy_all_levels.tsv` — all levels combined (renamed columns:
  `3-end mismatches`, `hits`, `missed`; renamed database labels: `mlas-MCRA`,
  `mlas-MCRA_CUR`, `MCRA_FULL`)
- `results/taxonomy/plots/sweep_taxonomy_all_levels.pdf` — combined multi-panel PDF:
  rows = taxonomic levels (phylum → genus), columns = databases + shared legend;
  y-axis = hits (%), x-axis = 3'-end allowed mismatches

```bash
Rscript scripts/parse_mismatch_sweep.R
```

All Step 6 scripts resolve paths relative to their own location and can be
called from any working directory.

**virtualPCR installation:** virtualPCR is available at
<https://github.com/rkalendar/virtualPCR>.  Clone the repository and build
the JAR once:

```bash
git clone https://github.com/rkalendar/virtualPCR.git
cd virtualPCR
mvn package -DskipTests          # requires Maven; produces dist/virtualPCR.jar
```

Then edit the `JAR=` line near the top of `run_virtualPCR.sh` to point to the
built JAR:

```bash
JAR="/path/to/virtualPCR/dist/virtualPCR.jar"
```

virtualPCR requires **OpenJDK ≥ 25**.  Install via conda (recommended):

```bash
conda create --name java25 openjdk=25 -c conda-forge --yes
```

The shell scripts activate this environment automatically.

---

## Dependencies

### R (≥ 4.3)

Install all required packages in one block:

```r
# CRAN packages
install.packages(c(
  "rentrez", "taxize", "stringr", "dplyr", "readr", "httr", "xml2",
  "ape", "phangorn", "ggplot2", "patchwork", "scales", "ggvenn",
  "pheatmap", "vegan", "pROC", "gridExtra", "tidyverse", "tidyr"
))

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Biostrings")
```

### External command-line tools

All tools are installable via conda (recommended):

```bash
conda install -c bioconda hmmer mafft blast cd-hit fasttree
conda install -c conda-forge ncbi-datasets-cli
```

| Tool | Version tested | Purpose |
|------|---------------|---------|
| MAFFT | 7.505 | Multiple sequence alignment |
| FastTree | 2.2.0 | Approximate ML tree (calibration step) |
| HMMER (hmmsearch) | 3.x | HMM-based gene verification and trimming |
| BLAST+ (blastn) | 2.x | Sequence search and annotation transfer |
| CD-HIT-EST | 4.x | Seed clustering |
| NCBI Datasets CLI | — | Genome batch download (optional) |

See `software_versions.tsv` for the exact versions used in this study.

### virtualPCR and Java (Step 6 only)

virtualPCR (<https://github.com/rkalendar/virtualPCR>) requires **OpenJDK ≥ 25**.
Clone and build the JAR (requires Maven), then point `run_virtualPCR.sh` to it
(see Step 6 above).  Install Java via conda (recommended):

```bash
conda create --name java25 openjdk=25 -c conda-forge --yes
```

The shell scripts activate this environment automatically.  If conda is not
available, any JDK ≥ 25 on `PATH` will work — remove or adapt the
`conda activate` block in `run_virtualPCR.sh`.

---

## Recommended workflow

```
Step 1  build_mcrA_db_ncbi_annot.R
          ↓ yang_2014_taxize.fasta
          ↓ mcrA_ncbi_genome_db.fasta

Step 2  calibrate_mcra_thresholds.R
          ↓ mcra_rank_thresholds.tsv

Step 3  build_mcrA_db.R
    (requires outputs of Steps 1 & 2)
          ↓ mcrA_ncbi_nt_db.fasta
          ↓ mcrA_ncbi_nt_cur_db.fasta

Step 4  compare_mcrA_databases2.R
    (requires outputs of Steps 1 & 3)
          ↓ comparison figures and tables

Step 5  convert_db_formats.R
    (requires outputs of Steps 1 & 3)
          ↓ <db>/dada2/, <db>/qiime2/, <db>/mothur/

Step 6  Primer coverage assessment (requires unzipped databases from Steps 1 & 3)

  # Unzip databases
  cd databases && for z in *.zip; do unzip "$z"; done && cd ..

  # Run mismatch sweep (0–3 mismatches)
  bash scripts/run_mismatch_sweep.sh

  # Parse sweep results → tables + multi-panel PDF
  Rscript scripts/parse_mismatch_sweep.R

  # Optional: parse single-run results (N3_ERRORS=1) separately
  bash scripts/run_virtualPCR.sh
  Rscript scripts/parse_virtualPCR.R
          ↓ results/primer_coverage/
          ↓ results/primer_coverage_sweep_errorsN/
          ↓ results/mismatch_sweep_summary.tsv
          ↓ results/taxonomy/sweep_taxonomy_all_levels.tsv
          ↓ results/taxonomy/plots/sweep_taxonomy_all_levels.pdf
```

Steps 1–5 read and write files in the current working directory; run each
script from the folder containing its input files.  Step 6 scripts resolve
all paths relative to their own location and can be called from any directory.

---

## Citation

If you use these databases or scripts, please cite:

- **This repository** (please use the GitHub URL / DOI once released)

- **Angel R, Claus P, Conrad R** (2012) Activation of methanogenesis in arid
  biological soil crusts despite the absence of liquid water. *ISME Journal*
  **6**, 2476–2483. https://doi.org/10.1038/ismej.2011.141

- **Yang S, Liebner S, Alawi M, Ebenhöh O, Rennert T** (2014) Taxonomic
  database and cut-off value for processing mcrA gene 454 pyrosequencing data
  by MOTHUR and QIIME. *Journal of Microbiological Methods* **103**, 3–5.

- **Kalendar R *et al.*** (2024) virtualPCR — a tool for in silico primer
  testing. *Frontiers in Bioinformatics* **4**, 1464197.
  https://doi.org/10.3389/fbinf.2024.1464197

---

## License

This repository is released under the [MIT License](LICENSE).
Sequence data retrieved from NCBI are subject to the
[NCBI data use policies](https://www.ncbi.nlm.nih.gov/home/about/policies/).
