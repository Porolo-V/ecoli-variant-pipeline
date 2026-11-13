#!/usr/bin/env bash
# E. coli REL606 variant calling: reads -> annotated variants
# Usage: bash scripts/run_all.sh        (conda env "varcall" must be active)
#        OUT=results_check bash scripts/run_all.sh   (write to another folder)
set -euo pipefail

# ---------- Parameters ----------
SAMPLE=SRR2584863
REF_ID=NC_012967.1
THREADS=${THREADS:-4}
OUT=${OUT:-results}
MIN_LEN=50            # fastp: drop reads shorter than this after trimming
CUT_Q=20              # fastp: sliding-window mean quality cutoff
MIN_QUAL=30           # bcftools: minimum variant QUAL
MIN_DP=10             # bcftools: minimum depth
MIN_MQ=40             # bcftools: minimum mapping quality (removes repeat artefacts)

cd "$(dirname "$0")/.."
RAW=data/raw
REF=data/ref/REL606.fa
GBK=data/ref/REL606.gbk
ENA=https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/003/${SAMPLE}
mkdir -p "$RAW" data/ref "$OUT"/{qc,trimmed,aligned,variants}

# ---------- 1. Download ----------
if [ ! -s "$REF" ]; then
  efetch -db nuccore -id "$REF_ID" -format fasta > "$REF"
fi
if [ ! -s "$GBK" ]; then
  efetch -db nuccore -id "$REF_ID" -format gbwithparts > "$GBK"
fi
for i in 1 2; do
  f=$RAW/${SAMPLE}_${i}.fastq.gz
  [ -s "$f" ] || aria2c -x 8 -s 8 -k 1M -c --disable-ipv6=true --max-tries=0 \
      --retry-wait=5 -d "$RAW" -o "$(basename "$f")" "${ENA}/$(basename "$f")"
done

# ---------- 2. QC before trimming ----------
fastqc --outdir "$OUT/qc" -t 2 "$RAW"/${SAMPLE}_1.fastq.gz "$RAW"/${SAMPLE}_2.fastq.gz
multiqc "$OUT/qc" -o "$OUT/qc" --force

# ---------- 3. Trimming ----------
fastp \
  --in1 "$RAW/${SAMPLE}_1.fastq.gz" --in2 "$RAW/${SAMPLE}_2.fastq.gz" \
  --out1 "$OUT/trimmed/R1.fastq.gz" --out2 "$OUT/trimmed/R2.fastq.gz" \
  --detect_adapter_for_pe \
  --cut_right --cut_window_size 4 --cut_mean_quality "$CUT_Q" \
  --qualified_quality_phred 20 --length_required "$MIN_LEN" \
  --thread "$THREADS" \
  --html "$OUT/trimmed/fastp.html" --json "$OUT/trimmed/fastp.json"

# ---------- 4. Reference index ----------
[ -s "$REF.bwt.2bit.64" ] || bwa-mem2 index "$REF"
[ -s "$REF.fai" ] || samtools faidx "$REF"

# ---------- 5. Alignment ----------
bwa-mem2 mem -t "$THREADS" \
  -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA" \
  "$REF" "$OUT/trimmed/R1.fastq.gz" "$OUT/trimmed/R2.fastq.gz" \
  | samtools sort -@ "$THREADS" -o "$OUT/aligned/sample.sorted.bam" -
samtools index "$OUT/aligned/sample.sorted.bam"
samtools flagstat "$OUT/aligned/sample.sorted.bam" > "$OUT/aligned/flagstat.txt"

# ---------- 6. Remove duplicates ----------
samtools sort -n -@ "$THREADS" "$OUT/aligned/sample.sorted.bam" \
  | samtools fixmate -m - - \
  | samtools sort -@ "$THREADS" - \
  | samtools markdup -r - "$OUT/aligned/sample.dedup.bam"
samtools index "$OUT/aligned/sample.dedup.bam"
samtools flagstat "$OUT/aligned/sample.dedup.bam" > "$OUT/aligned/flagstat_dedup.txt"
samtools coverage "$OUT/aligned/sample.dedup.bam" > "$OUT/aligned/coverage.tsv"

# ---------- 7. Variant calling ----------
bcftools mpileup -Ou -f "$REF" "$OUT/aligned/sample.dedup.bam" \
  | bcftools call -mv --ploidy 1 -Oz -o "$OUT/variants/raw.vcf.gz"
bcftools index "$OUT/variants/raw.vcf.gz"

# ---------- 8. Filtering ----------
bcftools filter -i "QUAL>${MIN_QUAL} && INFO/DP>${MIN_DP} && INFO/MQ>${MIN_MQ}" -Oz \
  -o "$OUT/variants/filtered.vcf.gz" "$OUT/variants/raw.vcf.gz"
bcftools index "$OUT/variants/filtered.vcf.gz"
echo "Variants before filter: $(bcftools view -H "$OUT/variants/raw.vcf.gz" | wc -l)"
echo "Variants after filter:  $(bcftools view -H "$OUT/variants/filtered.vcf.gz" | wc -l)"

# ---------- 9. Annotation (snpEff) ----------
mkdir -p snpeff_data/REL606
cp "$GBK" snpeff_data/REL606/genes.gbk
printf 'data.dir = %s/snpeff_data\nREL606.genome : REL606\n' "$(pwd)" > snpeff.config
[ -s snpeff_data/REL606/snpEffectPredictor.bin ] || \
  snpEff build -genbank -c snpeff.config -v REL606
snpEff -c snpeff.config -stats "$OUT/variants/snpEff_summary.html" REL606 \
  "$OUT/variants/filtered.vcf.gz" > "$OUT/variants/annotated.vcf"

# ---------- 10. Variant table ----------
bcftools query -f '%POS\t%REF\t%ALT\t%QUAL\t%INFO/DP\t%INFO/MQ\t%INFO/ANN\n' \
  "$OUT/variants/annotated.vcf" \
| awk -F'\t' 'BEGIN{OFS="\t"; print "POS","REF","ALT","QUAL","DP","MQ","EFFECT","IMPACT","GENE","PROTEIN_CHANGE","FLAG"}
  {
    split($7,e,","); split(e[1],a,"|");
    flag = ($5 < 30) ? "low_depth" : "";
    print $1,$2,$3,$4,$5,$6,a[2],a[3],a[4],a[11],flag
  }' > "$OUT/variants/variants_table.tsv"

echo "Done. Table: $OUT/variants/variants_table.tsv"
