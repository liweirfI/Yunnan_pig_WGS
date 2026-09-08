# Yunnan pig WGS — analysis workflow

Whole-genome sequencing workflow for *A whole-genome variant dataset for 234 indigenous pigs from
nine breeds in Yunnan, China*.

- Raw reads: NCBI SRA [PRJNA1510310](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1510310)
- Variants: EMBL-EBI EVA [PRJEB125594](https://www.ebi.ac.uk/eva/?eva-study=PRJEB125594)
- Reference: *Sus scrofa* Sscrofa11.1 (GCA_000003025.6 / GCF_000003025.6)

## Scripts

| Script | Step | Output |
| --- | --- | --- |
| `00_config.sh` | Paths, tool locations and thresholds — the only file to edit | — |
| `01_fastp_fastqc.sh` | Read filtering (fastp) and quality assessment (FastQC) | Table 1 |
| `02_mapping_bwa_mem2.sh` | `bwa-mem2 mem -M` piped into `samtools sort`, then read groups | BAM |
| `03_markdup_bqsr.sh` | Duplicate removal (Sambamba), optional BQSR | Analysis-ready BAM |
| `04_bam_qc_depth.sh` | `samtools flagstat` and `mosdepth` | Table 2 |
| `05_haplotypecaller_gvcf.sh` | GATK HaplotypeCaller in `-ERC GVCF` mode, 1 Mb windows | Per-sample GVCF |
| `06_combine_genotype.sh` | GATK CombineGVCFs and GenotypeGVCFs | Raw cohort VCF |
| `07_hardfilter_vcftools.sh` | GATK VariantFiltration, then VCFtools population filtering | Released SNP/InDel sets |
| `08_annovar_annotate.sh` | ANNOVAR against Ensembl release 104 | Figure 4 |
| `09_variant_stats.sh` | Ti/Tv, mutation spectrum, per-chromosome density | Figure 2, Figure 3, Table 3 |
| `10_report_actual_parameters.sh` | Collects the command lines each tool recorded in its own output | Provenance report |
| `run_all.sh` | Driver for all steps, or a subset | — |

`example/` holds a sample sheet (`samples.tsv`) and the Sscrofa11.1 chromosome list
(`chr_length.tsv`) in the formats the scripts expect.

## Requirements

fastp 0.23.4 · FastQC 0.12.1 · BWA-MEM2 2.2.1 · SAMtools 1.17 · Sambamba 1.0.1 · mosdepth 0.3.6 ·
GATK 4.4.0.0 · BCFtools 1.17 · VCFtools 0.1.16 · ANNOVAR 2019Oct24 (with UCSC `gtfToGenePred`) ·
Python 3.6+ (standard library only)

## Usage

Point `00_config.sh` at your reference genome, sample sheet and tool paths, then:

```bash
bash run_all.sh          # every step
bash run_all.sh 01 02 03 # selected steps
```

Every step detects finished output and skips it, so an interrupted run can be restarted by
re-issuing the same command. The published dataset was produced by running these same commands
with samples scheduled in parallel across a compute cluster.

## Read filtering parameters

`01_fastp_fastqc.sh` invokes fastp with `-w 7 -z 7` and otherwise relies on the defaults of
fastp v0.23.4, which are:

- adapters trimmed by per-pair overlap analysis (paired-end reads need no adapter sequence list)
- `-q 15 -u 40` — discard a read when more than 40% of its bases fall below Phred quality 15
- `-n 5` — discard a read containing more than 5 ambiguous (N) bases
- `-l 15` — discard reads shorter than 15 bp after trimming

Changing `FASTP_OPTS` produces a different clean read set from the one deposited under
PRJNA1510310. `10_report_actual_parameters.sh` reads the command line fastp records in each of its
JSON reports, so the settings behind the deposited data can always be confirmed directly.

## Base quality score recalibration

GATK BQSR requires a known-sites VCF. No curated truth set of variants exists for these breeds,
which is why the variant set is hard-filtered instead. `03_markdup_bqsr.sh` runs BaseRecalibrator
and ApplyBQSR when `KNOWN_SITES` is set in `00_config.sh`, and otherwise passes the
duplicate-removed BAM straight to variant calling.

## Citation

See `CITATION.cff`. Please cite both this workflow and the associated data descriptor.

## License

MIT — see `LICENSE`.
