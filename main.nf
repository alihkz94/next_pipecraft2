#!/usr/bin/env nextflow

nextflow.enable.dsl=2


params.outdir  = "$projectDir/results"
forward_ch= Channel.of(params.forward_primer)
reverse_ch= Channel.of(params.reverse_primer)


log.info """\
    PIPECRAFT2
    ===================================
    reads        : ${params.reads}
    outdir       : ${params.outdir}
    """
    .stripIndent()


process FASTQC {
    label 'FASTQC'
    tag "FASTQC on $sample_id"
    publishDir params.outdir

    input:
    tuple val(sample_id), path(reads)

    output:
    path "fastqc_${sample_id}_logs"

    script:
    """
    mkdir fastqc_${sample_id}_logs
    fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads}
    """
}

process MultiQC {
    label 'MULTIQC'
    publishDir params.outdir, mode:'copy' 

    input:
    path '*'

    output:
    path 'multiqc_report.html'

    script:
    """
    multiqc .
    """
}

process convertPrimer{
    publishDir params.outdir, mode:'copy'

    input:
    val forward_ch
    val reverse_ch
    
    output:
    path '**.fasta'

    shell:
    """
    Primer.sh ${params.forward_primer} ${params.reverse_primer}
    """
} 


workflow {
    Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .set { read_pairs_ch }

    fastqc_ch = FASTQC(read_pairs_ch)
 /*   MULTIQC(quant_ch.mix(fastqc_ch).collect()) */
    fastqc_ch
            .collect()
            .flatten()
            .view()
            .set {single_ch}
    MultiQC(single_ch)
    primer_ch= convertPrimer(forward_ch,reverse_ch) | collect | flatten | view
}





