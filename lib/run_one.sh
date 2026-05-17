#!/bin/bash
# ═══════════════════════════════════════════════════════════════
# hermes-h2/lib/run_one.sh — 单类型三软件遗传力
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

mkdir -p "$OUTDIR"/{filtered,plink,gcta,gemma,ldak}
LOG="$OUTDIR/${LABEL}.log"
exec > >(tee -a "$LOG") 2>&1

echo "══════════════════════════════════════════════"
echo "  hermes-h2 — $LABEL  |  $(date)"
echo "  MAF>$MAF  missing<$GENO  trait_col=$TRAIT_COL"
echo "══════════════════════════════════════════════"

# ══ 1. Filter VCF ══
echo "[1/6] Filtering VCF..."
FILT_VCF="$OUTDIR/filtered/${LABEL}.filtered.vcf.gz"
$BCFTOOLS view -i "F_MISSING < $GENO" -q ${MAF}:minor \
    -Oz -o "$FILT_VCF" "$VCF" --threads 4
$BCFTOOLS index -f "$FILT_VCF"

# ══ 1.5. Split multiallelic ══
echo "[1.5/6] Splitting multiallelic variants..."
$BCFTOOLS norm -m-any "$FILT_VCF" | \
    $BCFTOOLS view -m2 -M2 -Oz -o "${FILT_VCF/.vcf.gz/.biallelic.vcf.gz}"
$BCFTOOLS index -f "${FILT_VCF/.vcf.gz/.biallelic.vcf.gz}"
FILT_VCF="${FILT_VCF/.vcf.gz/.biallelic.vcf.gz}"

# ══ 2. VCF → PLINK ══
echo "[2/6] Converting to PLINK..."
$PLINK --vcf "$FILT_VCF" --set-all-var-ids @:# \
    --make-bed --out "$OUTDIR/plink/${LABEL}_tmp" --threads 4 --silent

# Sort samples to match phenotype order
python3 - "$PHENO_CSV" "$OUTDIR" << 'PYEOF'
import sys, csv
pheno_csv, outdir = sys.argv[1], sys.argv[2]
with open(pheno_csv) as f:
    reader = csv.reader(f)
    header = next(reader)
    ids = [row[0] for row in reader]
with open(f"{outdir}/sample_order.txt", 'w') as f:
    for iid in ids:
        f.write(f"0 {iid}\n")
PYEOF

$PLINK --bfile "$OUTDIR/plink/${LABEL}_tmp" \
    --keep "$OUTDIR/sample_order.txt" \
    --indiv-sort f "$OUTDIR/sample_order.txt" \
    --make-bed --out "$OUTDIR/plink/${LABEL}" --threads 4 --silent

rm -f "$OUTDIR"/plink/${LABEL}_tmp.* "$OUTDIR"/sample_order.txt
N_VAR=$(wc -l < "$OUTDIR/plink/${LABEL}.bim")
N_SAM=$(wc -l < "$OUTDIR/plink/${LABEL}.fam")
echo "  $N_VAR variants × $N_SAM samples"

# ══ 3. Phenotype ══
echo "[3/6] Preparing phenotype..."
python3 - "$PHENO_CSV" "$OUTDIR" "$LABEL" "$TRAIT_COL" << 'PYEOF'
import sys, csv
pheno_csv, outdir, label, trait_col = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])

fam_ids = []
with open(f"{outdir}/plink/{label}.fam") as f:
    for l in f:
        fam_ids.append(l.strip().split()[1])

with open(pheno_csv) as f:
    reader = csv.reader(f)
    header = next(reader)
    phe_dict = {}
    for row in reader:
        if len(row) > trait_col and row[trait_col] != 'NA':
            phe_dict[row[0]] = row[trait_col + 1]

# GCTA/LDAK: FID=0, IID=id, trait, header, missing=-9
with open(f"{outdir}/gcta/pheno.txt", 'w') as f:
    f.write(f"FID IID {label}\n")
    for iid in fam_ids:
        f.write(f"0 {iid} {phe_dict.get(iid, '-9')}\n")

# GEMMA: single column, matched to FAM order
with open(f"{outdir}/gemma/pheno_single.txt", 'w') as f:
    for iid in fam_ids:
        f.write(f"{phe_dict.get(iid, 'NA')}\n")

n_valid = sum(1 for v in phe_dict.values())
print(f"  {n_valid} phenotyped, {len(fam_ids)} total")
PYEOF

# ══ 4. GCTA ══
echo "[4/6] GCTA REML..."
$GCTA --bfile "$OUTDIR/plink/${LABEL}" --make-grm-bin \
    --out "$OUTDIR/gcta/grm" --threads 8 2>/dev/null || true
$GCTA --reml --grm "$OUTDIR/gcta/grm" --pheno "$OUTDIR/gcta/pheno.txt" \
    --out "$OUTDIR/gcta/h2" --threads 8 2>/dev/null || true

if [ -f "$OUTDIR/gcta/h2.hsq" ]; then
    GCTA_H2=$(grep "V(G)/Vp" "$OUTDIR/gcta/h2.hsq" | awk '{printf "%.4f", $2}')
    GCTA_SE=$(grep "V(G)/Vp" "$OUTDIR/gcta/h2.hsq" | awk '{printf "%.4f", $3}')
    echo "  GCTA: $GCTA_H2 ± $GCTA_SE"
else
    GCTA_H2="NA"; GCTA_SE="NA"; echo "  GCTA: FAILED"
fi

# ══ 5. GEMMA ══
echo "[5/6] GEMMA REML..."
awk '{$6=0; print}' OFS='\t' "$OUTDIR/plink/${LABEL}.fam" > "$OUTDIR/gemma/${LABEL}.fam"
cp "$OUTDIR/plink/${LABEL}.bed" "$OUTDIR/gemma/${LABEL}.bed"
cp "$OUTDIR/plink/${LABEL}.bim" "$OUTDIR/gemma/${LABEL}.bim"

cd "$OUTDIR/gemma"
$GEMMA -bfile "${LABEL}" -gk 2 -o "${LABEL}_kin" 2>/dev/null || true

GEMMA_OUT=$($GEMMA -bfile "${LABEL}" -k "output/${LABEL}_kin.sXX.txt" \
    -lmm 1 -miss 1.0 -maf 0 -p pheno_single.txt -o "${LABEL}_h2" 2>/dev/null) || true

GEMMA_H2=$(echo "$GEMMA_OUT" | grep "pve estimate" | awk -F'=' '{gsub(/ /,""); printf "%.4f", $2}')
GEMMA_SE=$(echo "$GEMMA_OUT" | grep "se(pve)" | awk -F'=' '{gsub(/ /,""); printf "%.4f", $2}')
if [ -n "$GEMMA_H2" ]; then
    echo "  GEMMA: $GEMMA_H2 ± $GEMMA_SE"
else
    GEMMA_H2="NA"; GEMMA_SE="NA"; echo "  GEMMA: FAILED"
fi
cd - > /dev/null

# ══ 6. LDAK ══
echo "[6/6] LDAK REML..."
cp "$OUTDIR/plink/${LABEL}.bed" "$OUTDIR/ldak/${LABEL}.bed"
cp "$OUTDIR/plink/${LABEL}.fam" "$OUTDIR/ldak/${LABEL}.fam"

if [ "$TYPE" != "SNP" ]; then
    awk 'BEGIN{FS=OFS="\t"} {if(NF>=6){a1=substr($5,1,1);a2=substr($6,1,1);
        if(a1==a2){a1="A";a2="T"};$5=a1;$6=a2};print}' \
        "$OUTDIR/plink/${LABEL}.bim" > "$OUTDIR/ldak/${LABEL}.bim"
else
    cp "$OUTDIR/plink/${LABEL}.bim" "$OUTDIR/ldak/${LABEL}.bim"
fi

cd "$OUTDIR/ldak"
$LDAK --calc-kins-direct "${LABEL}_kin" --bfile "${LABEL}" --power -0.25 2>/dev/null || true
$LDAK --reml "${LABEL}_h2" --grm "${LABEL}_kin" --pheno ../gcta/pheno.txt 2>/dev/null || true

if [ -f "${LABEL}_h2.reml" ]; then
    LDAK_H2=$(grep "Her_K1" "${LABEL}_h2.reml" | awk '{printf "%.4f", $2}')
    LDAK_SE=$(grep "Her_K1" "${LABEL}_h2.reml" | awk '{printf "%.4f", $3}')
    echo "  LDAK: $LDAK_H2 ± $LDAK_SE"
else
    LDAK_H2="NA"; LDAK_SE="NA"; echo "  LDAK: FAILED"
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
echo "Done: $(date)"
