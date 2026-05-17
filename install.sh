#!/bin/bash
# ═══════════════════════════════════════════════════════════════
# hermes-h2/install.sh — 自动下载 GCTA/GEMMA/LDAK
# ═══════════════════════════════════════════════════════════════
set -euo pipefail
BINDIR="$(cd "$(dirname "$0")" && pwd)/bin"
mkdir -p "$BINDIR"
cd "$BINDIR"

echo "═══ Installing hermes-h2 dependencies ═══"

# ── GCTA v1.94.1 ──
if [ ! -f gcta64 ]; then
    echo "  Downloading GCTA..."
    curl -sLO "https://yanglab.westlake.edu.cn/software/gcta/bin/gcta-1.94.1-linux-kernel-4-x86_64.zip"
    unzip -qo gcta-1.94.1-linux-kernel-4-x86_64.zip
    cp gcta-1.94.1-linux-kernel-4-x86_64/gcta64 . 2>/dev/null || \
    cp gcta-1.94.1-linux-kernel-4-x86_64/gcta-1.94.1 . 2>/dev/null
    chmod +x gcta64
    rm -rf gcta-1.94.1* __MACOSX
    echo "  GCTA: $(./gcta64 2>&1 | head -1)"
fi

# ── GEMMA v0.98.5 ──
if [ ! -f gemma ]; then
    echo "  Downloading GEMMA..."
    curl -sLO "https://github.com/genetics-statistics/GEMMA/releases/download/v0.98.5/gemma-0.98.5-linux-static-AMD64.gz"
    gunzip -f gemma-0.98.5-linux-static-AMD64.gz
    mv gemma-0.98.5-linux-static-AMD64 gemma
    chmod +x gemma
    echo "  GEMMA: $(./gemma 2>&1 | head -1)"
fi

# ── LDAK v6.2 ──
if [ ! -f ldak6.2 ]; then
    echo "  Downloading LDAK..."
    curl -sLO "https://raw.githubusercontent.com/dougspeed/LDAK/main/ldak6.2.linux"
    mv ldak6.2.linux ldak6.2
    chmod +x ldak6.2
    echo "  LDAK: $(./ldak6.2 2>&1 | grep "Version" | head -1)"
fi

# ── bcftools + plink2 via conda (optional) ──
echo ""
echo "  Note: bcftools and plink2 are also required."
echo "  Install via: conda install -c bioconda bcftools plink2"
echo ""
echo "═══ Done ═══"
