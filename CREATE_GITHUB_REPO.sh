#!/bin/bash

# MS2 Metabolite-Color Phenotype Analysis Pipeline
# GitHub Repository Creation Script

echo "╔════════════════════════════════════════════════════════════════╗"
echo "║  Creating GitHub Repository: Rhodotorula_MS2_pheno_explore    ║"
echo "╚════════════════════════════════════════════════════════════════╝"
echo ""

# Check if gh is installed
if ! command -v gh &> /dev/null; then
    echo "❌ ERROR: GitHub CLI (gh) is not installed"
    echo ""
    echo "Install it with:"
    echo "  macOS:   brew install gh"
    echo "  Linux:   https://github.com/cli/cli/releases"
    echo "  Windows: choco install gh"
    echo ""
    echo "After installation, run this script again."
    exit 1
fi

# Check if authenticated
if ! gh auth status &> /dev/null; then
    echo "❌ ERROR: Not authenticated with GitHub"
    echo ""
    echo "Run: gh auth login"
    echo "Then run this script again."
    exit 1
fi

echo "✓ GitHub CLI ready"
echo "✓ Authenticated as: $(gh auth status | head -1)"
echo ""

# Create the repository
echo "Creating repository: stajichlab/Rhodotorula_MS2_pheno_explore"
gh repo create stajichlab/Rhodotorula_MS2_pheno_explore \
  --public \
  --source=. \
  --remote=origin \
  --push \
  --description "MS2 metabolite-color phenotype analysis pipeline for Rhodotorula strains"

if [ $? -eq 0 ]; then
    echo ""
    echo "╔════════════════════════════════════════════════════════════════╗"
    echo "║                    ✅ REPOSITORY CREATED!                      ║"
    echo "╚════════════════════════════════════════════════════════════════╝"
    echo ""
    echo "Repository URL: https://github.com/stajichlab/Rhodotorula_MS2_pheno_explore"
    echo ""
    echo "Next steps:"
    echo "  1. Visit the repository on GitHub"
    echo "  2. Add topics: bioinformatics, metabolomics, yeast, rhodotorula, color-phenotype"
    echo "  3. Enable GitHub Discussions"
    echo "  4. Set up branch protection on 'main'"
    echo "  5. Add collaborators as needed"
    echo ""
    echo "Documentation:"
    echo "  • Start: analysis/00_START_HERE.md"
    echo "  • Findings: analysis/FINDINGS_SUMMARY.md"
    echo "  • Methods: analysis/README.md"
    echo ""
else
    echo ""
    echo "❌ Repository creation failed"
    echo ""
    echo "Troubleshooting:"
    echo "  • Check that 'stajichlab' organization exists and you have permission"
    echo "  • Or use: gh repo create Rhodotorula_MS2_pheno_explore --public"
    echo "  • Then manually transfer to stajichlab organization"
    exit 1
fi
