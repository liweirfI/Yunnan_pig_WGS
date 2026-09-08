#!/usr/bin/env bash
# 05_haplotypecaller_gvcf.sh -- single-sample variant calling in GVCF mode.
# Methods: "Post-alignment Processing and Variant Calling".
# Each sample is called in 1 Mb windows (-L) so the 234 samples finish in wall-clock time;
# the per-sample windows are then concatenated back into one GVCF per sample.
source "$(dirname "$0")/00_config.sh"

CHUNK=${CHUNK:-1000000}

read_samples | while IFS=$'\t' read -r sample fq_dir r1 r2 breed; do
    gvcf=$D_GVCF/${sample}.g.vcf.gz
    [[ -s $gvcf ]] && continue
    chunk_dir=$D_GVCF/chunks/$sample
    mkdir -p "$chunk_dir"
    : > "$chunk_dir/list.txt"

    while IFS=$'\t' read -r chr len alias order; do
        for start in $(seq 1 "$CHUNK" "$len"); do
            end=$(( start + CHUNK - 1 )); (( end > len )) && end=$len
            out=$chunk_dir/${chr}_${start}_${end}.g.vcf.gz
            echo "$out" >> "$chunk_dir/list.txt"
            [[ -s $out ]] && continue
            "$GATK" --java-options "-Xmx10G" HaplotypeCaller \
                -R "$REF" -I "$D_CLEANBAM/${sample}.analysis_ready.bam" \
                -L "${chr}:${start}-${end}" \
                --read-filter GoodCigarReadFilter \
                --native-pair-hmm-threads 8 --sample-ploidy 2 \
                --emit-ref-confidence GVCF \
                -O "$out"
        done
    done < <(grep -v '^#' "$CHRLIST" | awk 'NF>=2')

    "$BCFTOOLS" concat --threads "$THREADS" -f "$chunk_dir/list.txt" -Oz -o "$gvcf"
    "$BCFTOOLS" index -t --threads "$THREADS" "$gvcf"
done

echo "[05] done -> $D_GVCF/<sample>.g.vcf.gz"
