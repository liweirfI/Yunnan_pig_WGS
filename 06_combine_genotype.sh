#!/usr/bin/env bash
# 06_combine_genotype.sh -- CombineGVCFs across the 234 samples, then joint genotyping.
# Methods: "Post-alignment Processing and Variant Calling".
# Both steps run per chromosome; GenotypeGVCFs output is concatenated into the raw cohort VCF
# (~62,228,159 raw variants).
source "$(dirname "$0")/00_config.sh"

gvcf_args=()
while IFS=$'\t' read -r sample fq_dir r1 r2 breed; do
    gvcf_args+=( -V "$D_GVCF/${sample}.g.vcf.gz" )
done < <(read_samples)

: > "$D_VCF/chr_vcf.list"
while IFS=$'\t' read -r chr len alias order; do
    combined=$D_COMBINE/cohort_${chr}.g.vcf.gz
    raw=$D_VCF/raw_${chr}.vcf.gz
    echo "$raw" >> "$D_VCF/chr_vcf.list"

    [[ -s $combined ]] || "$GATK" --java-options "-Xmx${JAVA_MEM}" CombineGVCFs \
        -R "$REF" -L "$chr" "${gvcf_args[@]}" -O "$combined"

    [[ -s $raw ]] || "$GATK" --java-options "-Xmx${JAVA_MEM}" GenotypeGVCFs \
        -R "$REF" -L "$chr" -V "$combined" \
        --tmp-dir "$TMP_DIR" -O "$raw"
done < <(grep -v '^#' "$CHRLIST" | awk 'NF>=2')

RAW_VCF=$D_VCF/yunnan234_raw.vcf.gz
if [[ ! -s $RAW_VCF ]]; then
    "$BCFTOOLS" concat --threads "$THREADS" -f "$D_VCF/chr_vcf.list" -Oz -o "$RAW_VCF"
    "$BCFTOOLS" index -t --threads "$THREADS" "$RAW_VCF"
fi

echo "raw variants: $("$BCFTOOLS" index -n "$RAW_VCF")"
echo "[06] done -> $RAW_VCF"
