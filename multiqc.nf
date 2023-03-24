#!/usr/bin/env nextflow

nextflow.enable.dsl=2


params.projectDir = '.'
params.reads   = "$projectDir/data/F_{1,2}.fastq"
params.multiqc = "$projectDir/multiqc"
params.outdir  = "$projectDir/results"



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
}



//channel.fromPath('data/*.{fq, fastq.gz}', hidden: true)
//    .set {reads_ch}

// process fastqc {
//    input:
//    path reads
//
//    script:
//    """
//    echo $reads
//    """
//}  


//workflow{
//    def reads =channel.fromPath('data/*.{fq, fastq.gz}')
//    fastqc(reads)
//}




process fastqc {
    input:
    path reads

    script:
    """
    cat -n 2 $reads
    """
}

workflow {
    def reads = channel.fromPath('data/*.fq').view()
    fastqc(reads)
}


/* PRIMER TRIMMING */

primer_forward = Channel.of ()





