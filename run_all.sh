#!/usr/bin/env bash
# run_all.sh — full replication pipeline for "Drive-Level Momentum in the NFL"
#
# Reproduces every figure, table, and headline number in momentum.tex
# starting either from the live nflfastR feed (default) or from the cached
# drives.csv (skip step 1 if drives.rds is already present).
#
# Usage: bash run_all.sh

set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$PROJECT_ROOT"
mkdir -p replication_logs

echo "=== Step 1: Build drive-level dataset ==="
Rscript code/01_build_drives.R 2>&1 | tee replication_logs/01_build_drives.log

echo "=== Step 2: Main analysis (regressions, figures, tables) ==="
Rscript code/02_analyze_momentum.R 2>&1 | tee replication_logs/02_analyze_momentum.log

echo "=== Step 3: Robustness gauntlet ==="
Rscript code/03_robustness.R 2>&1 | tee replication_logs/03_robustness.log

echo "=== Step 4: Capture sessionInfo ==="
R --vanilla --quiet -e 'sessionInfo()' > replication_logs/sessionInfo.txt 2>&1

echo "=== Step 5: Compile paper_v2.tex ==="
if command -v pdflatex >/dev/null 2>&1; then
  pdflatex -interaction=nonstopmode paper_v2.tex > replication_logs/pdflatex_pass1.log 2>&1
  bibtex paper_v2 > replication_logs/bibtex.log 2>&1 || true
  pdflatex -interaction=nonstopmode paper_v2.tex > replication_logs/pdflatex_pass2.log 2>&1
  pdflatex -interaction=nonstopmode paper_v2.tex > replication_logs/pdflatex_pass3.log 2>&1
  echo "Compiled: paper_v2.pdf"
else
  echo "pdflatex not found; skipping LaTeX compile."
fi

echo "PIPELINE_OK" > replication_logs/STATUS
echo "Done. See replication_logs/ for run outputs."
