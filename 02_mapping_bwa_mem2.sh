#!/usr/bin/env bash
# 02_mapping_bwa_mem2.sh -- align clean reads to Sscrofa11.1 and add read groups.
# Methods: "Data Preprocessing and Mapping".  BWA-MEM2 mem -M piped straight into samtools sort,
# so no intermediate SAM file is written.
source "$(dirname "$0")/00_config.sh"

[[ -f ${REF}.bwt.2bit.64 ]] || "$BWA_MEM2" index "$REF"
[[ -f ${REF}.fai ]]         || "$SAMTOOLS" faidx "$REF"
[[ -f ${REF%.*}.dict ]]     || "$GATK" CreateSequenceDictionary -R "$REF"

read_samples | while IFS=$'\t' read -r sample fq_dir r1 r2 breed; do
    bam=$D_BAM/${sample}.sort.bam
    rgbam=$D_BAM/${sample}.sort.addRG.bam
    [[ -s $rgbam ]] && continue

    "$BWA_MEM2" mem -t "$THREADS" -M "$REF" \
        "$D_FASTP/${sample}_qc.R1.fq.gz" "$D_FASTP/${sample}_qc.R2.fq.gz" \
      | "$SAMTOOLS" sort -@ "$THREADS" -T "$TMP_DIR/$sample" -O bam - > "$bam"

    # ID / LB / PL / SM / PU as described in the Methods.  One library was built per animal, so
    # LB carries the sample ID; set it to the batch or BioProject if that is not true for you.
    "$SAMTOOLS" addreplacerg -@ "$THREADS" \
        -r "ID:${sample}" -r "LB:${sample}" -r "PL:DNBSEQ" \
        -r "SM:${sample}" -r "PU:dnbseq_t7" \
        -o "$rgbam" "$bam"
    "$SAMTOOLS" index -@ "$THREADS" "$rgbam"
    rm -f "$bam"
done

echo "[02] done -> $D_BAM/<sample>.sort.addRG.bam"
