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

workflow {
    reads_ch = channel.fromFilePairs(params.reads, checkIfExists: true)
    FASTP(reads_ch)
}
