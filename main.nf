process FASTP {
    publishDir "${params.outdir}/trimmed", mode: 'copy'

    input:
    tuple val(id), path(reads)

    output:
    tuple val(id), path("${id}_R1.trim.fastq.gz"), path("${id}_R2.trim.fastq.gz"), emit: reads
    path "${id}.fastp.json", emit: json
    path "${id}.fastp.html"

    script:
    """
    fastp \\
      --in1 ${reads[0]} --in2 ${reads[1]} \\
      --out1 ${id}_R1.trim.fastq.gz --out2 ${id}_R2.trim.fastq.gz \\
      --detect_adapter_for_pe \\
      --cut_right --cut_window_size 4 --cut_mean_quality 20 \\
      --qualified_quality_phred 20 --length_required 50 \\
      --thread ${task.cpus} \\
      --json ${id}.fastp.json --html ${id}.fastp.html
    """
}

process BWA_INDEX {
    input:
    path ref

    output:
    tuple path(ref), path("${ref}.*"), emit: index

    script:
    """
    bwa-mem2 index ${ref}
    """
}

process ALIGN {
    publishDir "${params.outdir}/aligned", mode: 'copy'

    input:
    tuple val(id), path(r1), path(r2)
    tuple path(ref), path(index_files)

    output:
    tuple val(id), path("${id}.sorted.bam"), path("${id}.sorted.bam.bai"), emit: bam
    path "${id}.flagstat.txt"

    script:
    """
    bwa-mem2 mem -t ${task.cpus} \\
      -R '@RG\\tID:${id}\\tSM:${id}\\tPL:ILLUMINA' \\
      ${ref} ${r1} ${r2} \\
      | samtools sort -@ ${task.cpus} -o ${id}.sorted.bam -
    samtools index ${id}.sorted.bam
    samtools flagstat ${id}.sorted.bam > ${id}.flagstat.txt
    """
}
process MARKDUP {
    publishDir "${params.outdir}/aligned", mode: 'copy'

    input:
    tuple val(id), path(bam), path(bai)

    output:
    tuple val(id), path("${id}.dedup.bam"), path("${id}.dedup.bam.bai"), emit: bam
    path "${id}.dedup.flagstat.txt"
    path "${id}.coverage.tsv"

    script:
    """
    samtools sort -n -@ ${task.cpus} ${bam} \\
      | samtools fixmate -m - - \\
      | samtools sort -@ ${task.cpus} - \\
      | samtools markdup -r - ${id}.dedup.bam
    samtools index ${id}.dedup.bam
    samtools flagstat ${id}.dedup.bam > ${id}.dedup.flagstat.txt
    samtools coverage ${id}.dedup.bam > ${id}.coverage.tsv
    """
}

process CALL_VARIANTS {
    publishDir "${params.outdir}/variants", mode: 'copy'

    input:
    tuple val(id), path(bam), path(bai)
    tuple path(ref), path(index_files)
    path fai

    output:
    tuple val(id), path("${id}.filtered.vcf.gz"), path("${id}.filtered.vcf.gz.csi"), emit: vcf
    path "${id}.raw.vcf.gz"

    script:
    """
    bcftools mpileup -Ou -f ${ref} ${bam} \\
      | bcftools call -mv --ploidy 1 -Oz -o ${id}.raw.vcf.gz
    bcftools index ${id}.raw.vcf.gz
    bcftools filter \\
      -i 'QUAL>${params.min_qual} && INFO/DP>${params.min_dp} && INFO/MQ>${params.min_mq}' \\
      -Oz -o ${id}.filtered.vcf.gz ${id}.raw.vcf.gz
    bcftools index ${id}.filtered.vcf.gz
    echo "raw: \$(bcftools view -H ${id}.raw.vcf.gz | wc -l)  filtered: \$(bcftools view -H ${id}.filtered.vcf.gz | wc -l)"
    """
}

process ANNOTATE {
    publishDir "${params.outdir}/variants", mode: 'copy'

    input:
    tuple val(id), path(vcf), path(csi)
    path gbk

    output:
    path "${id}.annotated.vcf"
    path "${id}.variants_table.tsv"
    path "${id}.snpEff_summary.html"

    script:
    """
    mkdir -p snpeff_data/REL606
    cp ${gbk} snpeff_data/REL606/genes.gbk
    printf 'data.dir = %s/snpeff_data\\nREL606.genome : REL606\\n' "\$(pwd)" > snpeff.config
    snpEff build -genbank -c snpeff.config REL606

    snpEff -c snpeff.config -stats ${id}.snpEff_summary.html REL606 ${vcf} > ${id}.annotated.vcf

    bcftools query -f '%POS\\t%REF\\t%ALT\\t%QUAL\\t%INFO/DP\\t%INFO/MQ\\t%INFO/ANN\\n' ${id}.annotated.vcf \\
    | awk -F'\\t' 'BEGIN{OFS="\\t"; print "POS","REF","ALT","QUAL","DP","MQ","EFFECT","IMPACT","GENE","PROTEIN_CHANGE","FLAG"}
      {
        split(\$7,e,","); split(e[1],a,"|");
        flag = (\$5 < 30) ? "low_depth" : "";
        print \$1,\$2,\$3,\$4,\$5,\$6,a[2],a[3],a[4],a[11],flag
      }' > ${id}.variants_table.tsv
    """
}

workflow {
    reads_ch = channel.fromFilePairs(params.reads, checkIfExists: true)
    ref_ch   = channel.fromPath(params.ref, checkIfExists: true)
    fai_ch   = channel.fromPath("${params.ref}.fai", checkIfExists: true)
    gbk_ch   = channel.fromPath(params.gbk, checkIfExists: true)

    FASTP(reads_ch)
    BWA_INDEX(ref_ch)
    ALIGN(FASTP.out.reads, BWA_INDEX.out.index)
    MARKDUP(ALIGN.out.bam)
    CALL_VARIANTS(MARKDUP.out.bam, BWA_INDEX.out.index, fai_ch)
    ANNOTATE(CALL_VARIANTS.out.vcf, gbk_ch)
}
