nextflow.enable.dsl = 2
params.threads = params.threads ?: 15
params.species_name = params.species_name ?: 'unknown_species'

process NANOCHECK1 {
    publishDir "results/nanoplot_before_trim", mode: 'copy'
    input:
    path fastq_file
    output:
    path 'NanoPlot_CHECK_1', emit: nanoplot_out
    script:
    """
    NanoPlot -t ${params.threads} --fastq $fastq_file --tsv_stats -o NanoPlot_CHECK_1
    """
}

process TRIM {
    publishDir "results/trimmed_fastq", mode: 'copy'
    input:
    path fastq_file
    output:
    path "${params.species_name}.trimmed.fastq", emit: trimmed_fastq
    script:
    """
    NanoFilt -q 10 $fastq_file > ${params.species_name}.trimmed.fastq
    """
}

process NANOCHECK2 {
    publishDir "results/nanoplot_after_trim", mode: 'copy'
    input:
    path trimmed_fastq
    output:
    path 'NanoPlot_CHECK_2', emit: nanoplot_out2
    script:
    """
    NanoPlot -t ${params.threads} --fastq $trimmed_fastq --tsv_stats -o NanoPlot_CHECK_2
    """
}

workflow {
    Channel.fromPath(params.fastfiles).set { fastq_ch }

    fastq_ch.view { "Input FASTQ: $it" }

    NANOCHECK1(fastq_ch)
    trimmed = TRIM(fastq_ch)
    NANOCHECK2(trimmed.trimmed_fastq)
}
