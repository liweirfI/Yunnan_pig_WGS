#!/usr/bin/env bash
# 03_markdup_bqsr.sh -- remove PCR duplicates (Sambamba), optionally recalibrate base qualities
# (GATK BQSR), and produce the analysis-ready BAM.
# Methods: "Post-alignment Processing and Variant Calling".
# BQSR runs only when KNOWN_SITES is set in 00_config.sh -- GATK needs a truth VCF of known
# variants, and there is no curated one for these breeds.  See README note 2.
source "$(dirname "$0")/00_config.sh"

read_samples | while IFS=$'\t' read -r sample fq_dir r1 r2 breed; do
    rgbam=$D_BAM/${sample}.sort.addRG.bam
    dedup=$D_CLEANBAM/${sample}.rmDup.bam
    final=$D_CLEANBAM/${sample}.analysis_ready.bam
    [[ -s $final ]] && continue

    if [[ ! -s $dedup ]]; then
        "$SAMBAMBA" markdup -r -t "$THREADS" --tmpdir="$TMP_DIR" \
            --overflow-list-size 1000000 --hash-table-size 1000000 \
            "$rgbam" "$dedup"
    fi

    if [[ -n $KNOWN_SITES ]]; then
        "$GATK" --java-options "-Xmx${JAVA_MEM}" BaseRecalibrator \
            -R "$REF" -I "$dedup" --known-sites "$KNOWN_SITES" \
            -O "$D_CLEANBAM/${sample}.recal.table"
        "$GATK" --java-options "-Xmx${JAVA_MEM}" ApplyBQSR \
            -R "$REF" -I "$dedup" \
            --bqsr-recal-file "$D_CLEANBAM/${sample}.recal.table" \
            -O "$final"
    else
        ln -sf "$dedup" "$final"
        "$SAMTOOLS" index -@ "$THREADS" "$final"
    fi
done

echo "[03] done -> $D_CLEANBAM/<sample>.analysis_ready.bam"
