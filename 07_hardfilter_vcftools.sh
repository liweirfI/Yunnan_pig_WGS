#!/usr/bin/env bash
# 07_hardfilter_vcftools.sh -- GATK hard filtering followed by population-level filtering.
# Methods: "Variant Filtration and Annotation".
# Yields the released call set: biallelic SNPs and InDels with MAF > 0.01 and missing rate < 0.1.
source "$(dirname "$0")/00_config.sh"

RAW_VCF=$D_VCF/yunnan234_raw.vcf.gz
MAX_MISSING=$("$PYTHON" -c "print(1 - $MAX_MISS)")   # vcftools counts the kept fraction

for TYPE in SNP INDEL; do
    sel=$D_FILTER/raw_${TYPE}.vcf.gz
    flg=$D_FILTER/hardfilter_${TYPE}.vcf.gz
    [[ $TYPE == SNP ]] && EXPR=$SNP_FILTER || EXPR=$INDEL_FILTER

    [[ -s $sel ]] || "$GATK" --java-options "-Xmx${JAVA_MEM}" SelectVariants \
        -R "$REF" -V "$RAW_VCF" --select-type-to-include "$TYPE" -O "$sel"

    # VariantFiltration only flags; --remove-filtered-all below drops the flagged records.
    [[ -s $flg ]] || "$GATK" --java-options "-Xmx${JAVA_MEM}" VariantFiltration \
        -R "$REF" -V "$sel" \
        --filter-name "${TYPE}_hard_filter" --filter-expression "$EXPR" \
        -O "$flg"

    # Population-level filtering (VCFtools): biallelic only, MAF and genotype missingness.
    final=$D_FILTER/yunnan234_${TYPE}.vcf.gz
    if [[ ! -s $final ]]; then
      "$VCFTOOLS" --gzvcf "$flg" \
        --remove-filtered-all \
        --min-alleles 2 --max-alleles 2 \
        --maf "$MAF" --max-missing "$MAX_MISSING" \
        --recode --recode-INFO-all \
        --stdout | "$BCFTOOLS" view -Oz -o "$final"
      "$BCFTOOLS" index -t --threads "$THREADS" "$final"
    fi

    echo "${TYPE}: $("$BCFTOOLS" index -n "$D_FILTER/yunnan234_${TYPE}.vcf.gz") retained"
done

# Combined release VCF (the file deposited at EVA).
if [[ ! -s $D_FILTER/yunnan234_SNP_INDEL.vcf.gz ]]; then
"$BCFTOOLS" concat --threads "$THREADS" -a -Oz \
    -o "$D_FILTER/yunnan234_SNP_INDEL.vcf.gz" \
    "$D_FILTER/yunnan234_SNP.vcf.gz" "$D_FILTER/yunnan234_INDEL.vcf.gz"
"$BCFTOOLS" index -t --threads "$THREADS" "$D_FILTER/yunnan234_SNP_INDEL.vcf.gz"
fi

echo "[07] done -> $D_FILTER/yunnan234_{SNP,INDEL,SNP_INDEL}.vcf.gz"
