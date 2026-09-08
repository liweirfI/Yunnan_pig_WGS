#!/usr/bin/env bash
# 09_variant_stats.sh -- Ti/Tv ratio, SNP mutation spectrum and per-chromosome variant density.
# Technical Validation: Figure 2 (mutation types, Ti/Tv = 2.46) and Figure 3 / Table 3
# (SNP and InDel density along the genome).
source "$(dirname "$0")/00_config.sh"

for TYPE in SNP INDEL; do
    tsv=$D_STATS/${TYPE}_chrom_pos.tsv
    [[ -s $tsv ]] || "$BCFTOOLS" query -f '%CHROM\t%POS\t%REF\t%ALT\n' \
        "$D_FILTER/yunnan234_${TYPE}.vcf.gz" > "$tsv"
done

"$PYTHON" - "$D_STATS" "$CHRLIST" <<'PY'
import collections, os, sys
stats_dir, chrlist = sys.argv[1], sys.argv[2]

# ---- Figure 2: six strand-collapsed mutation classes + Ti/Tv ----
comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
spectrum = collections.Counter()
snp_per_chr = collections.Counter()
with open(os.path.join(stats_dir, "SNP_chrom_pos.tsv")) as fh:
    for line in fh:
        chrom, pos, ref, alt = line.rstrip("\n").split("\t")
        snp_per_chr[chrom] += 1
        if ref in comp and alt in comp:
            if ref in ("A", "G"):          # collapse to a pyrimidine-first reference
                ref, alt = comp[ref], comp[alt]
            spectrum[f"{ref}>{alt}"] += 1
ti = spectrum["C>T"] + spectrum["T>C"]
tv = sum(v for k, v in spectrum.items() if k not in ("C>T", "T>C"))
total = sum(spectrum.values())
if not total:
    sys.exit("no biallelic SNP records found - did step 07 run?")
titv = f"{ti/tv:.4f}" if tv else "NA"     # tv == 0 only on a degenerate call set
with open(os.path.join(stats_dir, "Fig2_mutation_spectrum.tsv"), "w") as fh:
    fh.write("Mutation\tCount\tPercent\n")
    for k, v in sorted(spectrum.items(), key=lambda x: -x[1]):
        fh.write(f"{k}\t{v}\t{100*v/total:.2f}\n")
    fh.write(f"Ti/Tv\t{titv}\t\n")
print(f"SNPs {total} | Ti {ti} | Tv {tv} | Ti/Tv {titv}")

# ---- Figure 3 / Table 3: per-chromosome counts and density ----
indel_per_chr = collections.Counter()
with open(os.path.join(stats_dir, "INDEL_chrom_pos.tsv")) as fh:
    for line in fh:
        indel_per_chr[line.split("\t", 1)[0]] += 1
sizes = {}
order = []
for line in open(chrlist):
    f = line.split()
    if len(f) >= 2 and not line.startswith("#"):
        sizes[f[0]] = int(f[1]); order.append(f[0])
with open(os.path.join(stats_dir, "Table3_variant_density.tsv"), "w") as fh:
    fh.write("Chromosome\tSize (Mb)\tSNP number\tSNP density (bp/SNP)\t"
             "InDel number\tInDel density (bp/InDel)\n")
    for c in order:
        n_s, n_i, L = snp_per_chr[c], indel_per_chr[c], sizes[c]
        d_s = f"{L/n_s:.1f}" if n_s else "NA"   # a chromosome may carry no called variant
        d_i = f"{L/n_i:.1f}" if n_i else "NA"   # (e.g. Y, or a chromosome left out of the run)
        fh.write(f"{c}\t{L/1e6:.2f}\t{n_s}\t{d_s}\t{n_i}\t{d_i}\n")
print("wrote Fig2_mutation_spectrum.tsv and Table3_variant_density.tsv")
PY

echo "[09] done -> $D_STATS/"
