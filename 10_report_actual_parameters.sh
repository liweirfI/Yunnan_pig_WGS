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

section "fastp (version and command line recorded in the JSON reports)"
# The JSON records the version that actually processed the reads, which is the number the
# methods section needs; the installed binary may since have been upgraded.
"$PYTHON" - "$D_FASTP" <<'PY' >> "$REPORT"
import glob, json, os, sys
seen = set()
for js in sorted(glob.glob(os.path.join(sys.argv[1], "*.json"))):
    try:
        d = json.load(open(js))
        seen.add((d.get("summary", {}).get("fastp_version", "<not recorded>"),
                  d.get("command", "<no command field>")))
    except Exception as e:
        seen.add(("<unreadable>", f"{os.path.basename(js)}: {e}"))
if not seen:
    print("<no fastp JSON found>")
for version, cmd in sorted(seen):
    print(f"fastp {version}\n  {cmd}")
PY
echo "fastp binary installed now: $("$FASTP" --version 2>&1)" >> "$REPORT"

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

section "software versions of the binaries installed now (not necessarily those used above)"
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
