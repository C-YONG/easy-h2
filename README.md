# easy-h2

**Three-tool heritability estimation for SNP, INDEL, and SV.**

Estimates heritability using GCTA, GEMMA, and LDAK — three complementary methods — from filtered VCF and phenotype data.

## Install

```bash
git clone https://github.com/yourname/easy-h2.git
cd easy-h2
bash install.sh          # Downloads GCTA, GEMMA, LDAK
conda install -c bioconda bcftools plink2  # Prerequisites
```

## Quick start

```bash
# Single variant type
./easy-h2 single snp.vcf.gz pheno.csv SNP --out results/

# All three types
./easy-h2 batch snp.vcf.gz indel.vcf.gz sv.vcf.gz pheno.csv --out results/
```

## Input

| File | Format |
|------|--------|
| VCF | bgzip-compressed, indexed |
| Phenotype CSV | `FID,trait1,trait2,...` with header, `NA` for missing |

## Output

```
results/
├── filtered/        Filtered VCFs
├── plink/           PLINK binary files
├── gcta/            GCTA REML output
├── gemma/           GEMMA REML output
├── ldak/            LDAK REML output
└── summary.tsv      Final table
```

`summary.tsv`:

| Type | GCTA_h2 | GCTA_SE | GEMMA_h2 | GEMMA_SE | LDAK_h2 | LDAK_SE |
|------|---------|---------|----------|----------|---------|---------|
| SNP  | 0.319   | 0.093   | 0.320    | 0.104    | 0.343   | 0.095   |
| INDEL| 0.356   | 0.100   | 0.350    | 0.112    | 0.388   | 0.102   |
| SV   | 0.308   | 0.084   | 0.256    | 0.081    | 0.296   | 0.086   |

## Options

```
--out DIR       Output directory (default: ./easy-h2-out)
--trait N       Trait column (0-indexed after FID, default: 1)
```

## Default filtering

| Type  | MAF  | Missingness |
|-------|------|-------------|
| SNP   | 0.05 | 0.05        |
| INDEL | 0.05 | 0.05        |
| SV    | 0.01 | 0.20        |

## Dependencies

Auto-downloaded by `install.sh`:
- GCTA v1.94.1
- GEMMA v0.98.5
- LDAK v6.2

Required via conda/apt:
- bcftools ≥1.10
- plink2 ≥2.0
- python3

## Cite

If you use easy-h2, please cite the underlying tools:
- Yang et al. (2011) GCTA. *AJHG*
- Zhou & Stephens (2012) GEMMA. *Nature Genetics*
- Speed et al. (2017) LDAK. *Nature Genetics*
