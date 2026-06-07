#!/bin/bash
# =============================================================================
# LithoInv -- Publish to GitHub and EarthArXiv
# =============================================================================
# Run this script after reviewing the initial commit:
#   bash publish.sh
# =============================================================================

set -e

echo "========================================"
echo " LithoInv -- Publication Script"
echo "========================================"

# --- GitHub ---
echo ""
echo "[1/3] GitHub repository setup"
echo "----------------------------------------"
echo "Manual steps (do these in your browser first):"
echo "  1. Go to https://github.com/new"
echo "  2. Repository name: LithoInv"
echo "  3. Description: Lithology-constrained 3D joint inversion of gravity and magnetics"
echo "  4. Set to PUBLIC"
echo "  5. Do NOT initialize with README (we already have one)"
echo "  6. Click 'Create repository'"
echo ""
echo "After creating the repo, run these commands:"
echo ""
echo "  git remote add origin https://github.com/YOUR_USERNAME/LithoInv.git"
echo "  git push -u origin master"
echo ""
echo "Replace YOUR_USERNAME with your GitHub username."

# --- EarthArXiv ---
echo ""
echo "[2/3] EarthArXiv preprint upload"
echo "----------------------------------------"
echo "Manual steps:"
echo "  1. Go to https://eartharxiv.org/"
echo "  2. Create an account or log in"
echo "  3. Click 'Submit a Preprint'"
echo "  4. Select category: Earth Sciences -- Geophysics"
echo "  5. Upload preprint.tex (and any figures if added later)"
echo "  6. Fill in metadata:"
echo "     Title: LithoInv: Lithology-constrained joint gravity and magnetic"
echo "            inversion via GMM clustering and boundary optimization"
echo "     Author: Thomas Bagley"
echo "     Keywords: geophysics, inversion, gravity, magnetics, lithology, open source"
echo "     License: MIT"
echo "  7. Submit. You will receive a DOI after acceptance."
echo ""
echo "  After acceptance, update the DOI in:"
echo "    - README.md (How to Cite section)"
echo "    - preprint.tex (if recompiling)"
echo "    - runner.cpp (Citation line in banner)"

# --- Verify ---
echo ""
echo "[3/3] Verification checklist"
echo "----------------------------------------"
echo "  [ ] LICENSE file present (MIT, Copyright 2026 Thomas Bagley)"
echo "  [ ] NOTICE file present"
echo "  [ ] README.md updated with description, citation, build instructions"
echo "  [ ] preprint.tex compiled and uploaded"
echo "  [ ] Git repo pushed to GitHub (public)"
echo "  [ ] EarthArXiv preprint submitted"
echo "  [ ] Banner prints from InversionRunner constructor (runner.cpp)"
echo ""
echo "After EarthArXiv acceptance, update all placeholders:"
echo "  README.md:   https://eartharxiv.org/repository/view/XXXXX/"
echo "  preprint.tex: https://github.com/XXXXX/LithoInv"
echo "========================================"
