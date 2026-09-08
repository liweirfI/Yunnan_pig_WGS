#!/usr/bin/env bash
# 08_annovar_annotate.sh -- functional annotation with ANNOVAR against Ensembl release 104.
# Methods: "Variant Filtration and Annotation".  Produces the category counts behind Figure 4.
# ANNOVAR ships no pig database, so a refGene-style database is built from the Ensembl GTF.
source "$(dirname "$0")/00_config.sh"

BUILD=${BUILD:-susScr11}
DB=$D_ANNO/db
mkdir -p "$DB"

# ---- one-off: GTF -> GenePred -> ANNOVAR refGene database ----
if [[ ! -s $DB/${BUILD}_refGeneMrna.fa ]]; then
    "$ANNOVAR_GTF2GP" -genePredExt "$GTF" "$DB/${BUILD}_refGene0.txt"
    nl -ba "$DB/${BUILD}_refGene0.txt" > "$DB/${BUILD}_refGene.txt"
    perl "$(dirname "$TABLE_ANNOVAR")/retrieve_seq_from_fasta.pl" \
        "$DB/${BUILD}_refGene.txt" --seqfile "$REF" --format refGene \
        -out "$DB/${BUILD}_refGeneMrna.fa"
fi

for TYPE in SNP INDEL; do
    vcf=$D_FILTER/yunnan234_${TYPE}.vcf.gz
    prefix=$D_ANNO/yunnan234_${TYPE}
    [[ -s ${prefix}.${BUILD}_multianno.txt ]] || \
        perl "$TABLE_ANNOVAR" "$vcf" "$DB" \
            -buildver "$BUILD" -out "$prefix" \
            -protocol refGene -operation g \
            -vcfinput -remove -polish -nastring . -thread "$THREADS"
done

# ---- Figure 4: counts and percentages per genomic category ----
"$PYTHON" - "$D_ANNO" "$BUILD" "$D_STATS/Fig4_annotation_categories.tsv" <<'PY'
import collections, sys, os
anno_dir, build, out = sys.argv[1], sys.argv[2], sys.argv[3]
table = {}
for t in ("SNP", "INDEL"):
    path = os.path.join(anno_dir, f"yunnan234_{t}.{build}_multianno.txt")
    counts = collections.Counter()
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        col = header.index("Func.refGene")
        for line in fh:
            counts[line.split("\t")[col]] += 1
    table[t] = counts
cats = sorted(set(table["SNP"]) | set(table["INDEL"]),
              key=lambda c: -table["SNP"][c])
with open(out, "w") as fh:
    fh.write("Category\tSNP_count\tSNP_percent\tInDel_count\tInDel_percent\n")
    ts, ti = sum(table["SNP"].values()), sum(table["INDEL"].values())
    for c in cats:
        fh.write(f"{c}\t{table['SNP'][c]}\t{100*table['SNP'][c]/ts:.2f}\t"
                 f"{table['INDEL'][c]}\t{100*table['INDEL'][c]/ti:.2f}\n")
    print(f"SNP {ts} / InDel {ti} annotated")
PY

echo "[08] done -> $D_STATS/Fig4_annotation_categories.tsv"
