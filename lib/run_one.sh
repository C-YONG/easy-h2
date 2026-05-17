#!/bin/bash
# ═══════════════════════════════════════════════════════════════
# easy-h2/lib/run_one.sh — 单类型三软件遗传力
# ═══════════════════════════════════════════════════════════════
set -euo pipefail

VCF="$1"; PHENO_CSV="$2"; TYPE="$3"; OUTDIR="$4"; TRAIT_COL="$5"
GCTA="$6"; GEMMA="$7"; LDAK="$8"
MAF_SNP="$9"; GENO_SNP="${10}"; MAF_INDEL="${11}"; GENO_INDEL="${12}"; MAF_SV="${13}"; GENO_SV="${14}"

case "$TYPE" in
    SNP|snp)   MAF=$MAF_SNP; GENO=$GENO_SNP; LABEL="SNP" ;;
    INDEL|indel) MAF=$MAF_INDEL; GENO=$GENO_INDEL; LABEL="INDEL" ;;
    SV|sv)     MAF=$MAF_SV; GENO=$GENO_SV; LABEL="SV" ;;
    *) echo "ERROR: type must be SNP/INDEL/SV"; exit 1 ;;
esac

# ── Per-type subdirectory (no overlap in batch mode) ──
TDIR="$OUTDIR/$LABEL"
mkdir -p "$TDIR"/{filtered,plink,gcta,gemma,ldak}
LOG="$TDIR/${LABEL}.log"
exec > >(tee -a "$LOG") 2>&1

echo "══════════════════════════════════════════════"
echo "  easy-h2 — $LABEL  |  $(date)"
echo "  MAF>$MAF  missing<$GENO  trait_col=$TRAIT_COL"
echo "══════════════════════════════════════════════"

ERRORS=0

# ══ 1. Filter VCF ══
echo "[1/6] Filtering VCF..."
FILT_VCF="$TDIR/filtered/${LABEL}.filtered.vcf.gz"
$BCFTOOLS view -i "F_MISSING < $GENO" -q ${MAF}:minor -m2 -M2 \
    -Oz -o "$FILT_VCF" "$VCF" --threads 4
$BCFTOOLS index -f "$FILT_VCF"

# ══ 2. VCF → PLINK ══
echo "[2/6] Converting to PLINK..."
# Variant IDs: SNP use chr:pos:ref:alt, INDEL/SV use chr:pos (alleles too long)
if [ "$TYPE" = "SNP" ]; then
    ID_FMT="@:#:\$r:\$a"
else
    ID_FMT="@:#"
fi
$PLINK --vcf "$FILT_VCF" --set-all-var-ids $ID_FMT \
    --allow-extra-chr \
    --make-bed --out "$TDIR/plink/${LABEL}_tmp" --threads 4 --silent

# Sort samples to match phenotype order
python3 - "$PHENO_CSV" "$TDIR" << 'PYEOF'
import sys, csv
pheno_csv, tdir = sys.argv[1], sys.argv[2]
with open(pheno_csv) as f:
    reader = csv.reader(f)
    next(reader)
    ids = [row[0] for row in reader]
with open(f"{tdir}/sample_order.txt", 'w') as f:
    for iid in ids:
        f.write(f"0 {iid}\n")
PYEOF

$PLINK --bfile "$TDIR/plink/${LABEL}_tmp" \
    --keep "$TDIR/sample_order.txt" \
    --indiv-sort f "$TDIR/sample_order.txt" \
    --allow-extra-chr \
    --make-bed --out "$TDIR/plink/${LABEL}" --threads 4 --silent

rm -f "$TDIR"/plink/${LABEL}_tmp.* "$TDIR"/sample_order.txt
N_VAR=$(wc -l < "$TDIR/plink/${LABEL}.bim")
N_SAM=$(wc -l < "$TDIR/plink/${LABEL}.fam")
echo "  $N_VAR variants × $N_SAM samples"

# ══ 3. Phenotype ══
echo "[3/6] Preparing phenotype..."
python3 - "$PHENO_CSV" "$TDIR" "$LABEL" "$TRAIT_COL" << 'PYEOF'
import sys, csv
pheno_csv, tdir, label, trait_col = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])

fam_ids = []
with open(f"{tdir}/plink/{label}.fam") as f:
    for l in f:
        fam_ids.append(l.strip().split()[1])

with open(pheno_csv) as f:
    reader = csv.reader(f)
    next(reader)
    phe_dict = {}
    for row in reader:
        if len(row) > trait_col and row[trait_col] != 'NA':
            phe_dict[row[0]] = row[trait_col]

with open(f"{tdir}/gcta/pheno.txt", 'w') as f:
    f.write(f"FID IID {label}\n")
    for iid in fam_ids:
        f.write(f"0 {iid} {phe_dict.get(iid, '-9')}\n")

with open(f"{tdir}/gemma/pheno_single.txt", 'w') as f:
    for iid in fam_ids:
        f.write(f"{phe_dict.get(iid, 'NA')}\n")

n_valid = sum(1 for v in phe_dict.values())
print(f"  {n_valid} phenotyped, {len(fam_ids)} total")
PYEOF

# ══ 4. GCTA ══
echo "[4/6] GCTA REML..."
if $GCTA --bfile "$TDIR/plink/${LABEL}" --make-grm-bin \
    --out "$TDIR/gcta/grm" --threads 8 2>"$TDIR/gcta/gcta_grm.err" && \
   $GCTA --reml --grm "$TDIR/gcta/grm" --pheno "$TDIR/gcta/pheno.txt" \
    --out "$TDIR/gcta/h2" --threads 8 2>"$TDIR/gcta/gcta_reml.err"; then
    GCTA_H2=$(grep "V(G)/Vp" "$TDIR/gcta/h2.hsq" | awk '{printf "%.4f", $2}')
    GCTA_SE=$(grep "V(G)/Vp" "$TDIR/gcta/h2.hsq" | awk '{printf "%.4f", $3}')
    echo "  GCTA: $GCTA_H2 ± $GCTA_SE"
else
    GCTA_H2="NA"; GCTA_SE="NA"; ERRORS=$((ERRORS+1))
    echo "  GCTA: FAILED (see $TDIR/gcta/gcta_*.err)"
fi

# ══ 5. GEMMA ══
echo "[5/6] GEMMA REML..."
awk '{$6=0; print}' OFS='\t' "$TDIR/plink/${LABEL}.fam" > "$TDIR/gemma/${LABEL}.fam"
cp "$TDIR/plink/${LABEL}.bed" "$TDIR/gemma/${LABEL}.bed"
cp "$TDIR/plink/${LABEL}.bim" "$TDIR/gemma/${LABEL}.bim"

cd "$TDIR/gemma"
if $GEMMA -bfile "${LABEL}" -gk 2 -o "${LABEL}_kin" 2>gemma_kin.err; then
    GEMMA_OUT=$($GEMMA -bfile "${LABEL}" -k "output/${LABEL}_kin.sXX.txt" \
        -lmm 1 -miss 1.0 -maf 0 -p pheno_single.txt -o "${LABEL}_h2" 2>gemma_reml.err) || true
    GEMMA_H2=$(echo "$GEMMA_OUT" | grep "pve estimate" | awk -F'=' '{gsub(/ /,""); printf "%.4f", $2}')
    GEMMA_SE=$(echo "$GEMMA_OUT" | grep "se(pve)" | awk -F'=' '{gsub(/ /,""); printf "%.4f", $2}')
fi
if [ -n "${GEMMA_H2:-}" ]; then
    echo "  GEMMA: $GEMMA_H2 ± $GEMMA_SE"
else
    GEMMA_H2="NA"; GEMMA_SE="NA"; ERRORS=$((ERRORS+1))
    echo "  GEMMA: FAILED (see $TDIR/gemma/gemma_*.err)"
fi
cd - > /dev/null

# ══ 6. LDAK ══
echo "[6/6] LDAK REML..."
cp "$TDIR/plink/${LABEL}.bed" "$TDIR/ldak/${LABEL}.bed"
cp "$TDIR/plink/${LABEL}.fam" "$TDIR/ldak/${LABEL}.fam"

if [ "$TYPE" != "SNP" ]; then
    awk 'BEGIN{FS=OFS="\t"} {if(NF>=6){a1=substr($5,1,1);a2=substr($6,1,1);
        if(a1==a2){a1="A";a2="T"};$5=a1;$6=a2};print}' \
        "$TDIR/plink/${LABEL}.bim" > "$TDIR/ldak/${LABEL}.bim"
else
    cp "$TDIR/plink/${LABEL}.bim" "$TDIR/ldak/${LABEL}.bim"
fi

cd "$TDIR/ldak"
if $LDAK --calc-kins-direct "${LABEL}_kin" --bfile "${LABEL}" --power -0.25 2>ldak_kin.err && \
   $LDAK --reml "${LABEL}_h2" --grm "${LABEL}_kin" --pheno ../gcta/pheno.txt 2>ldak_reml.err; then
    LDAK_H2=$(grep "Her_K1" "${LABEL}_h2.reml" | awk '{printf "%.4f", $2}')
    LDAK_SE=$(grep "Her_K1" "${LABEL}_h2.reml" | awk '{printf "%.4f", $3}')
    echo "  LDAK: $LDAK_H2 ± $LDAK_SE"
else
    LDAK_H2="NA"; LDAK_SE="NA"; ERRORS=$((ERRORS+1))
    echo "  LDAK: FAILED (see $TDIR/ldak/ldak_*.err)"
fi
cd - > /dev/null

# ══ Summary ══
SUMMARY="$OUTDIR/summary.tsv"
if [ ! -f "$SUMMARY" ]; then
    printf "Type\tGCTA_h2\tGCTA_SE\tGEMMA_h2\tGEMMA_SE\tLDAK_h2\tLDAK_SE\n" > "$SUMMARY"
fi
printf "${LABEL}\t${GCTA_H2}\t${GCTA_SE}\t${GEMMA_H2}\t${GEMMA_SE}\t${LDAK_H2}\t${LDAK_SE}\n" >> "$SUMMARY"

echo ""
echo "═══ $LABEL: GCTA=$GCTA_H2 GEMMA=$GEMMA_H2 LDAK=$LDAK_H2 ═══"
if [ $ERRORS -gt 0 ]; then
    echo "WARNING: $ERRORS tool(s) failed — check *.err files"
fi
echo "Done: $(date)"
