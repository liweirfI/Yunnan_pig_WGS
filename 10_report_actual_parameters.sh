#!/usr/bin/env bash
# 10_report_actual_parameters.sh -- report the command lines that produced a given run.
# Run it on the machine that holds the output tree.  Every tool records its own invocation,
# so the parameters behind a dataset can be read back from the data itself:
#   fastp   -> "command" field of each JSON report
#   BWA-MEM2 / samtools / sambamba -> @PG lines in the BAM header
#   GATK    -> GATKCommandLine header lines in the GVCF and the raw VCF
#   VCFtools / ANNOVAR -> the log files they write next to their output
source "$(dirname "$0")/00_config.sh"
# This script only interrogates files that already exist.  A missing tool, an absent header
# line, or the SIGPIPE from "... | head -1" must not abort it, so strict mode is dropped here.
set +e +o pipefail

REPORT=$D_STATS/actual_parameters.txt
: > "$REPORT"

section() { printf '\n===== %s =====\n' "$1" >> "$REPORT"; }

section "fastp (unique command lines across all samples)"
"$PYTHON" - "$D_FASTP" <<'PY' >> "$REPORT"
import glob, json, os, sys
cmds = set()
for js in glob.glob(os.path.join(sys.argv[1], "*.json")):
    try:
        cmds.add(json.load(open(js)).get("command", "<no command field>"))
    except Exception as e:
        cmds.add(f"<unreadable {os.path.basename(js)}: {e}>")
print("\n".join(sorted(cmds)) or "<no fastp JSON found>")
PY
echo "fastp version: $("$FASTP" --version 2>&1)" >> "$REPORT"

section "BAM @PG (first sample)"
first_bam=$(ls "$D_CLEANBAM"/*.analysis_ready.bam 2>/dev/null | head -1)
[[ -n $first_bam ]] && "$SAMTOOLS" view -H "$first_bam" | grep '^@PG' >> "$REPORT"

section "GATK HaplotypeCaller"
first_gvcf=$(ls "$D_GVCF"/*.g.vcf.gz 2>/dev/null | head -1)
[[ -n $first_gvcf ]] && "$BCFTOOLS" view -h "$first_gvcf" | grep -m1 'GATKCommandLine' >> "$REPORT"

section "GATK GenotypeGVCFs / VariantFiltration"
[[ -s $D_VCF/yunnan234_raw.vcf.gz ]] && \
    "$BCFTOOLS" view -h "$D_VCF/yunnan234_raw.vcf.gz" | grep 'GATKCommandLine' >> "$REPORT"
[[ -s $D_FILTER/yunnan234_SNP.vcf.gz ]] && \
    "$BCFTOOLS" view -h "$D_FILTER/yunnan234_SNP.vcf.gz" | grep -E 'FILTER=<ID=|GATKCommandLine' >> "$REPORT"

section "software versions"
{
    "$BWA_MEM2" version 2>&1 | head -1
    "$SAMTOOLS" --version | head -1
    "$SAMBAMBA" --version 2>&1 | head -2
    "$GATK" --version 2>&1 | grep -i 'Genome Analysis Toolkit' | head -1
    "$BCFTOOLS" --version | head -1
    "$VCFTOOLS" --version 2>&1 | head -1
    "$FASTQC" --version 2>&1 | head -1
} >> "$REPORT"

cat "$REPORT"
echo "[10] done -> $REPORT"
exit 0
