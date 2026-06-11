#!/usr/bin/env bash
# run_virtualPCR.sh
#
# In silico PCR — mcrA primer coverage assessment using virtualPCR
#
# Primers:
#   Forward: mlas-mod-F  5'-GGYGGTGTMGGDTTCACMCARTA-3'
#   Reverse: mcrA-rev-R  5'-BGCGTAGTTVGGRTAGT-3'  (first 7 bases CGTTCAT excluded)
#
# Databases (relative to the repo root — unzip databases/*.zip first):
#   databases/mcrA_ncbi_genome_db/dada2/mcrA_ncbi_genome_db.fasta
#   databases/mcrA_ncbi_nt_cur_db/dada2/mcrA_ncbi_nt_cur_db.fasta
#   databases/mcrA_ncbi_nt_db/dada2/mcrA_ncbi_nt_db.fasta
#
# Requires: conda environment "java25" with OpenJDK >= 25
#   conda create --name java25 openjdk=25 -c conda-forge --yes
#
# Usage (from any directory):
#   bash scripts/run_virtualPCR.sh
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

# ── Paths (all relative to this script's location) ───────────────────────────
BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"   # scripts/
REPO_DIR="$(cd "${BASE_DIR}/.." && pwd)"                   # repo root
JAR="${BASE_DIR}/virtualPCR/dist/virtualPCR.jar"
RES_DIR="${REPO_DIR}/results/primer_coverage"
PRIMER_FILE="${RES_DIR}/primers.fasta"

# ── virtualPCR parameters ─────────────────────────────────────────────────────
N3_ERRORS=1
MIN_LEN=100
MAX_LEN=900

# ── Databases ─────────────────────────────────────────────────────────────────
# Unzip databases/*.zip into databases/ before running.
DATABASE_NAMES=(
  "mcrA_ncbi_genome_db"
  "mcrA_ncbi_nt_cur_db"
  "mcrA_ncbi_nt_db"
)
DATABASE_FASTAS=(
  "${REPO_DIR}/databases/mcrA_ncbi_genome_db/dada2/mcrA_ncbi_genome_db.fasta"
  "${REPO_DIR}/databases/mcrA_ncbi_nt_cur_db/dada2/mcrA_ncbi_nt_cur_db.fasta"
  "${REPO_DIR}/databases/mcrA_ncbi_nt_db/dada2/mcrA_ncbi_nt_db.fasta"
)

# ── Activate conda java25 ─────────────────────────────────────────────────────
CONDA_BASE="$(conda info --base 2>/dev/null)" || {
  echo "ERROR: conda not found on PATH. Activate your base conda environment first."
  exit 1
}
# shellcheck source=/dev/null
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate java25

echo "Java version: $(java -version 2>&1 | head -1)"
echo "JAR         : ${JAR}"
echo ""

# ── Write primer FASTA ────────────────────────────────────────────────────────
mkdir -p "${RES_DIR}"
cat > "${PRIMER_FILE}" <<'EOF'
>mlas-mod-F
GGYGGTGTMGGDTTCACMCARTA
>mcrA-rev-R
BGCGTAGTTVGGRTAGT
EOF
echo "Primer file : ${PRIMER_FILE}"
echo ""

# ── Run virtualPCR for each database ─────────────────────────────────────────
for i in "${!DATABASE_NAMES[@]}"; do
    DB_NAME="${DATABASE_NAMES[$i]}"
    FASTA="${DATABASE_FASTAS[$i]}"
    DB_OUT="${RES_DIR}/${DB_NAME}"
    CONFIG="${DB_OUT}/config.txt"
    REPORT="${DB_OUT}/report.out"

    echo "──────────────────────────────────────────────────────"
    echo "Database : ${DB_NAME}"
    echo "FASTA    : ${FASTA}"

    if [[ ! -f "${FASTA}" ]]; then
        echo "WARNING: FASTA not found — skipping."
        echo "  → Unzip databases/${DB_NAME}.zip first."
        continue
    fi

    mkdir -p "${DB_OUT}"
    [[ -f "${REPORT}" ]] && rm -f "${REPORT}"

    cat > "${CONFIG}" <<EOF
targets_path=${FASTA}
output_path=${REPORT}
primers_path=${PRIMER_FILE}
type=primer
linkedsearch=false
molecular=linear
number3errors=${N3_ERRORS}
minlen=${MIN_LEN}
maxlen=${MAX_LEN}
FRpairs=false
CTconversion=false
SequenceExtract=false
ShowPrimerAlignment=false
ShowOnlyAmplicons=false
ShowPCRProducts=true
ShowPrimerAlignmentPCRproduct=false
primerstatistic=false
EOF

    echo "Config   : ${CONFIG}"
    echo "Running virtualPCR..."
    java -jar "${JAR}" "${CONFIG}"
    echo "Done. Report: ${REPORT}"
    echo ""
done

echo "══════════════════════════════════════════════════════"
echo "All databases processed."
echo "Run: Rscript scripts/parse_virtualPCR.R"
echo "Results directory: ${RES_DIR}"
