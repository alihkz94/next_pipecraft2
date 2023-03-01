#!/usr/bin/env nextflow

params {
    exts = 'fasta,fasta.gz,fastq,fastq.gz,fa.gz,fa,fq,fq.gz,txt,txt.gz'
    min_seq_length = 1
    min_overlap = 1
    mismatch = 0
    no_indels = false
    seqs_to_keep = 'keep_all'
    discard_untrimmed = true
    output_dir = './output'
    cpus = 1
    mem = '1g'
}

nextflow.enable.dsl=2

process clipPrimers {
    input:
    file input_file from file_channel

    output:
    file("${output_dir}/${input_file.baseName}_trimmed.fasta") into trimmed_files_channel

    script:
    """
    # Remove primers from reads
    # Degenerate primers are allowed using IUPAC codes. Reverse complementary strings will also be searched.

    extension=${input_file.extension}
    mismatches="-e ${mismatch}"
    min_length="--minimum-length ${min_seq_length}"
    overlap="--overlap ${min_overlap}"
    cores="--cores ${cpus}"
    no_indels=${no_indels}
    discard_untrimmed=${discard_untrimmed}
    seqs_to_keep=${seqs_to_keep}

    # Prepare primers
    fwd_primer_array=( $(echo "${params.forward_primers}" | sed 's/,/ /g' | sed 's/I/N/g') )
    rev_primer_array=( $(echo "${params.reverse_primers}" | sed 's/,/ /g' | sed 's/I/N/g') )

    # Forward primer(s) to fasta file
    for ((i=1; i<=${#fwd_primer_array[@]}; i++)); do
        echo ">fwd_primer$i" >> tempdir2/fwd_primer.fasta
        echo ${fwd_primer_array[i-1]} >> tempdir2/fwd_primer.fasta
    done

    # Reverse primer(s) to fasta file
    for ((i=1; i<=${#rev_primer_array[@]}; i++)); do
        echo ">rev_primer$i" >> tempdir2/rev_primer.fasta
        echo ${rev_primer_array[i-1]} >> tempdir2/rev_primer.fasta
    done

    # Reverse complement REV primers
    seqkit seq --quiet -t dna -r -p tempdir2/rev_primer.fasta >> tempdir2/rev_primer_RC.fasta

    # Make linked primers files
    i=1
    while read LINE; do
        fwd_primer=$(echo $LINE | grep -v ">")
        if [ -z "$fwd_primer" ]; then
            :
        else
            while read LINE; do
                rev_primer_RC=$(echo $LINE | grep -v ">")
                if [ -z "$rev_primer_RC" ]; then
                    :
                else
                    echo ">primer$i" >> tempdir2/liked_fwd_revRC.fasta
                    echo "$fwd_primer...$rev_primer_RC" >> tempdir2/liked_fwd_revRC.fasta
                    ((i=i+1))
                fi
            done < tempdir2/rev_primer_RC.fasta
        fi
    done < tempdir2/fwd_primer.fasta

    # Clip primers with cutadapt
    if [[ $seqs_to_keep == "keep_all" ]]; then
        cutadapt --quiet --revcomp \
        $mismatches \
                $min_length \
        $overlap \
        ${no_indels:?'-n'} \
        $cores \
        ${discard_untrimmed:?'-u'} \
        -g file:tempdir2/liked_fwd_revRC.fasta \
        -g file:tempdir2/fwd_primer.fasta \
        -a file:tempdir2/rev_primer_RC.fasta \
        -o $output_dir/${input_file.baseName}_trimmed.fasta \
        $input_file

    elif [[ $seqs_to_keep == "keep_only_linked" ]]; then
        cutadapt --quiet --revcomp \
        $mismatches \
        $min_length \
        $overlap \
        ${no_indels:?'-n'} \
        $cores \
        ${discard_untrimmed:?'-u'} \
        -g file:tempdir2/liked_fwd_revRC.fasta \
        -o $output_dir/${input_file.baseName}_trimmed.fasta \
        $input_file
    fi

    # Clean up
    
    rm -f tempdir2/fwd_primer.fasta tempdir2/rev_primer.fasta tempdir2/rev_primer_RC.fasta tempdir2/liked_fwd_revRC.fasta
    """
}

channel.fromPath(params.exts)
    .map { file(it) }
    .set { file_channel }

workflow {
    trimmed_files_channel
        .filter { it != null }
        .into { trimmed_files }

    trimmed_files
        .collect()
        .set { trimmed_files_list }

    // print summary
    afterScript {
        println "\nExecution summary:"
        println " - Input files: ${file_channel.size()}"
        println " - Output files: ${trimmed_files_list.size()}"
    }
}



