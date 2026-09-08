#!/usr/bin/env bash
# 00_config.sh -- shared paths, tool versions and thresholds for the 234 Yunnan pig WGS dataset.
# Every step script starts with:  source "$(dirname "$0")/00_config.sh"
# Edit ONLY this file; the step scripts take no arguments.

set -euo pipefail

# ---------------- project layout ----------------
WORK=${WORK:-/disk192/users_dir/liw/2_3_yunnan_dir/1_2_yunnan234_dir}
SAMPLES=${SAMPLES:-$(dirname "${BASH_SOURCE[0]}")/example/samples.tsv}   # sample<TAB>fq_dir<TAB>R1<TAB>R2<TAB>breed
CHRLIST=${CHRLIST:-$(dirname "${BASH_SOURCE[0]}")/example/chr_length.tsv} # chr<TAB>length<TAB>alias<TAB>order

# ---------------- reference (Sscrofa11.1, GCF_000003025.6) ----------------
REF=${REF:-/disk192/users_dir/liw/reference_genomes_dir/pig/Sscrofa11.1.fa}
GTF=${GTF:-/disk192/users_dir/liw/reference_genomes_dir/pig/Sus_scrofa.Sscrofa11.1.104.gtf}
# Known-sites VCF for BQSR. Empty string = BQSR is skipped (see README, note 2).
KNOWN_SITES=${KNOWN_SITES:-}

# ---------------- tools (fill in absolute paths on your cluster) ----------------
FASTP=${FASTP:-fastp}                 # v0.23.4
FASTQC=${FASTQC:-fastqc}              # v0.12.1
BWA_MEM2=${BWA_MEM2:-bwa-mem2}        # v2.2.1
SAMTOOLS=${SAMTOOLS:-samtools}        # v1.17
SAMBAMBA=${SAMBAMBA:-sambamba}        # v1.0.1
MOSDEPTH=${MOSDEPTH:-mosdepth}        # v0.3.6
GATK=${GATK:-gatk}                    # v4.4.0.0
BCFTOOLS=${BCFTOOLS:-bcftools}        # v1.17
VCFTOOLS=${VCFTOOLS:-vcftools}        # v0.1.16
TABLE_ANNOVAR=${TABLE_ANNOVAR:-/opt/annovar/table_annovar.pl}       # v2019Oct24
ANNOVAR_GTF2GP=${ANNOVAR_GTF2GP:-/opt/annovar/gtfToGenePred}
PYTHON=${PYTHON:-python3}

# ---------------- runtime ----------------
THREADS=${THREADS:-30}
JAVA_MEM=${JAVA_MEM:-30g}
# Deliberately NOT called TMPDIR: that name is normally already exported by the shell, so
# ${TMPDIR:-...} would silently keep /tmp and fill it with sort and markdup spill files.
TMP_DIR=${TMP_DIR:-$WORK/tmp}

# ---------------- variant filtering (Methods: Variant Filtration) ----------------
SNP_FILTER=${SNP_FILTER:-"QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 3.0"}
INDEL_FILTER=${INDEL_FILTER:-"QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0"}
MAF=${MAF:-0.01}          # minor allele frequency threshold
MAX_MISS=${MAX_MISS:-0.1} # max genotype missing rate (vcftools --max-missing 0.9)

# ---------------- derived directories ----------------
D_FASTP=$WORK/1_fastp_dir
D_BAM=$WORK/2_bam_dir
D_CLEANBAM=$WORK/3_clean_bam_dir
D_QC=$WORK/3_bam_qc_dir
D_GVCF=$WORK/4_gVCF_dir
D_COMBINE=$WORK/5_combine_dir
D_VCF=$WORK/6_VCF_dir
D_FILTER=$WORK/7_filter_dir
D_ANNO=$WORK/8_annovar_dir
D_STATS=$WORK/9_stats_dir
mkdir -p "$TMP_DIR" "$D_FASTP" "$D_BAM" "$D_CLEANBAM" "$D_QC" "$D_GVCF" \
         "$D_COMBINE" "$D_VCF" "$D_FILTER" "$D_ANNO" "$D_STATS"

# Iterate the sample sheet: skips the header and blank lines.
read_samples() { grep -v '^#' "$SAMPLES" | awk 'NR>1 && NF>=5'; }
