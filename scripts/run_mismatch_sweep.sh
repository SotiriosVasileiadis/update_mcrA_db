#!/usr/bin/env bash
# run_mismatch_sweep.sh
#
# Runs run_virtualPCR.sh for each 3'-end mismatch level (0–3) by creating
# a temporary modified copy of the script — the original is never changed.
#
# Results land in:
#   results/primer_coverage_sweep_errors0/
#   results/primer_coverage_sweep_errors1/
#   results/primer_coverage_sweep_errors2/
#   results/primer_coverage_sweep_errors3/
#
# Usage (from any directory):
#   bash scripts/run_mismatch_sweep.sh
#   bash scripts/run_mismatch_sweep.sh 0 1 2    # run only specific levels
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"   # scripts/
REPO_DIR="$(cd "${BASE_DIR}/.." && pwd)"                   # repo root
ORIGINAL="${BASE_DIR}/run_virtualPCR.sh"

# Mismatch levels to test — override on the command line if desired
if [[ $# -eq 0 ]]; then
    ERRORS_LIST=(0 1 2 3)
else
    ERRORS_LIST=("$@")
fi

echo "══════════════════════════════════════════════════════"
echo "Mismatch sweep: number3errors = ${ERRORS_LIST[*]}"
echo "Original script: ${ORIGINAL}"
echo "══════════════════════════════════════════════════════"
echo ""

for NERR in "${ERRORS_LIST[@]}"; do

    RES_SUBDIR="results/primer_coverage_sweep_errors${NERR}"

    echo "──────────────────────────────────────────────────────"
    echo "Running with number3errors=${NERR}  →  ${RES_SUBDIR}/"
    echo "──────────────────────────────────────────────────────"

    TMP=$(mktemp "${TMPDIR:-/tmp}/run_vpcr_err${NERR}_XXXXXX")
    sed \
        -e "s|^N3_ERRORS=.*|N3_ERRORS=${NERR}|" \
        -e "s|^RES_DIR=.*|RES_DIR=\"\${REPO_DIR}/${RES_SUBDIR}\"|" \
        "${ORIGINAL}" > "${TMP}"
    chmod +x "${TMP}"

    echo "  Temp script : ${TMP}"
    bash "${TMP}"

    rm -f "${TMP}"
    echo "  Done → ${REPO_DIR}/${RES_SUBDIR}/"
    echo ""
done

echo "══════════════════════════════════════════════════════"
echo "Sweep complete."
echo "Run: Rscript scripts/parse_mismatch_sweep.R"
echo "══════════════════════════════════════════════════════"
