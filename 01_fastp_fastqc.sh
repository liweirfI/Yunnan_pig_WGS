#!/usr/bin/env bash
# 01_fastp_fastqc.sh -- raw-read filtering (fastp) and quality assessment (FastQC).
# Methods: "Data Preprocessing and Mapping".  Feeds Table 1 (raw/clean Gb, Q20, Q30, GC).
source "$(dirname "$0")/00_config.sh"

# The deposited dataset was filtered with fastp's defaults: only -w (worker threads) and
# -z (gzip level) were passed.  For v0.23.4 the defaults in force are:
#   adapters trimmed by per-pair overlap analysis (paired-end reads need no adapter list)
#   -q 15 -u 40   discard a read when >40% of its bases have Phred quality below 15
#   -n 5          discard a read containing more than 5 ambiguous (N) bases
#   -l 15         discard reads shorter than 15 bp after trimming
# Setting FASTP_OPTS yields a different clean read set from the one under PRJNA1510310.
FASTP_OPTS=${FASTP_OPTS:-}

read_samples | while IFS=$'\t' read -r sample fq_dir r1 r2 breed; do
    out1=$D_FASTP/${sample}_qc.R1.fq.gz
    out2=$D_FASTP/${sample}_qc.R2.fq.gz
    [[ -s $out2 ]] && continue
    "$FASTP" -w 7 -z 7 $FASTP_OPTS \
        -i "$fq_dir/${sample}${r1}" -I "$fq_dir/${sample}${r2}" \
        -o "$out1" -O "$out2" \
        -j "$D_FASTP/${sample}.json" -h "$D_FASTP/${sample}.html"
    "$FASTQC" --format fastq --threads 4 --outdir "$D_FASTP" "$out1" "$out2"
done

# ---- Table 1: per-sample raw/clean data volume, Q20, Q30, GC from the fastp JSON ----
"$PYTHON" - "$D_FASTP" "$D_STATS/Table1_sequencing_stats.tsv" <<'PY'
import json, sys, glob, os
fastp_dir, out = sys.argv[1], sys.argv[2]
rows = []
for js in sorted(glob.glob(os.path.join(fastp_dir, "*.json"))):
    d = json.load(open(js))
    before, after = d["summary"]["before_filtering"], d["summary"]["after_filtering"]
    rows.append([
        os.path.basename(js)[:-5],
        round(before["total_bases"] / 1e9, 4),          # Raw Data (G)
        round(after["total_bases"] / 1e9, 4),           # Clean Data (G)
        round(after["q20_rate"] * 100, 4),              # Q20 (%)
        round(after["q30_rate"] * 100, 4),              # Q30 (%)
        round(after["gc_content"] * 100, 2),            # GC (%)
    ])
with open(out, "w") as fh:
    fh.write("ID\tRaw Data (G)\tClean Data (G)\tQ20 (%)\tQ30 (%)\tGC Content (%)\n")
    for r in rows:
        fh.write("\t".join(map(str, r)) + "\n")
n = len(rows)
if n:
    m = [sum(r[i] for r in rows) / n for i in range(1, 6)]
    print(f"{n} samples | raw {m[0]:.2f} G | clean {m[1]:.2f} G | "
          f"Q20 {m[2]:.2f}% | Q30 {m[3]:.2f}% | GC {m[4]:.2f}%")
PY

echo "[01] done -> $D_STATS/Table1_sequencing_stats.tsv"
