#!/usr/bin/env bash
# run_all.sh -- run the whole pipeline end to end, or a subset of steps.
#   ./run_all.sh              run every step
#   ./run_all.sh 05 06 07     run only those steps
# Each step is restartable: finished outputs are detected and skipped.
set -euo pipefail
cd "$(dirname "$0")"

STEPS=(01_fastp_fastqc.sh
       02_mapping_bwa_mem2.sh
       03_markdup_bqsr.sh
       04_bam_qc_depth.sh
       05_haplotypecaller_gvcf.sh
       06_combine_genotype.sh
       07_hardfilter_vcftools.sh
       08_annovar_annotate.sh
       09_variant_stats.sh
       10_report_actual_parameters.sh)

if [[ $# -gt 0 ]]; then
    selected=()
    for want in "$@"; do
        for s in "${STEPS[@]}"; do [[ $s == ${want}_* ]] && selected+=("$s"); done
    done
    if [[ ${#selected[@]} -eq 0 ]]; then
        echo "no step matches: $*" >&2
        printf '  %s\n' "${STEPS[@]}" >&2
        exit 1
    fi
    STEPS=("${selected[@]}")
fi

mkdir -p logs
for s in "${STEPS[@]}"; do
    echo "=== $(date '+%F %T')  $s"
    bash "$s" 2>&1 | tee -a "logs/${s%.sh}.log"
done
echo "=== $(date '+%F %T')  finished: ${STEPS[*]}"
