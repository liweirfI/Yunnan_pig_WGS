#!/usr/bin/env bash
# 04_bam_qc_depth.sh -- per-sample sequencing depth, genome coverage and mapping rate.
# Technical Validation: feeds Table 2 (clean reads, mapped reads, mapping rate, 1x/4x coverage)
# and the reported mean depth of 11.67x (9.96x-17.76x).
source "$(dirname "$0")/00_config.sh"

read_samples | while IFS=$'\t' read -r sample fq_dir r1 r2 breed; do
    bam=$D_CLEANBAM/${sample}.analysis_ready.bam
    [[ -s $D_QC/${sample}.flagstat.txt ]] || \
        "$SAMTOOLS" flagstat -@ 10 "$bam" > "$D_QC/${sample}.flagstat.txt"
    [[ -s $D_QC/${sample}.mosdepth.summary.txt ]] || \
        "$MOSDEPTH" -t 5 -n --fast-mode "$D_QC/${sample}" "$bam"
done

# ---- Table 2 ----
"$PYTHON" - "$D_QC" "$D_STATS/Table2_mapping_stats.tsv" <<'PY'
import glob, gzip, os, re, sys
qc_dir, out = sys.argv[1], sys.argv[2]

def flagstat(path):
    total = mapped = 0
    for line in open(path):
        if " in total " in line and not total:
            total = int(line.split()[0])
        m = re.match(r"^(\d+) \+ \d+ mapped \(", line)
        if m:
            mapped = int(m.group(1))
    return total, mapped

def coverage(prefix):
    """1x and 4x genome coverage from the mosdepth global distribution."""
    cov = {}
    dist = prefix + ".mosdepth.global.dist.txt"
    if os.path.exists(dist):
        for line in open(dist):
            chrom, depth, frac = line.split()
            if chrom == "total":
                cov[int(float(depth))] = float(frac)
    mean = ""
    summ = prefix + ".mosdepth.summary.txt"
    if os.path.exists(summ):
        for line in open(summ):
            f = line.split()
            if f[0] == "total":
                mean = f[3]
    return mean, cov.get(1, ""), cov.get(4, "")

rows = []
for fs in sorted(glob.glob(os.path.join(qc_dir, "*.flagstat.txt"))):
    s = os.path.basename(fs)[:-len(".flagstat.txt")]
    total, mapped = flagstat(fs)
    mean, c1, c4 = coverage(os.path.join(qc_dir, s))
    rows.append([s, total, mapped,
                 round(100 * mapped / total, 4) if total else "",
                 mean,
                 round(100 * c1, 4) if c1 != "" else "",
                 round(100 * c4, 4) if c4 != "" else ""])
with open(out, "w") as fh:
    fh.write("ID\tClean reads\tMapped reads\tMapping rate (%)\tMean depth (x)\t"
             "Coverage (1x, %)\tCoverage (4x, %)\n")
    for r in rows:
        fh.write("\t".join(map(str, r)) + "\n")
if rows:
    mr = [r[3] for r in rows if r[3] != ""]
    dp = [float(r[4]) for r in rows if r[4] != ""]
    print(f"{len(rows)} samples | mapping rate {sum(mr)/len(mr):.2f}% | "
          f"depth mean {sum(dp)/len(dp):.2f}x range {min(dp):.2f}-{max(dp):.2f}x")
PY

echo "[04] done -> $D_STATS/Table2_mapping_stats.tsv"
