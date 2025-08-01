//Preprocessing: NanoPlot before/after trimming + NanoFilt.
nextflow.enable.dsl=2
workflow {
    Channel.fromPath(params.fastfiles, checkIfExists: true).set { fastq_ch }
NANOCHECK1(fastq_ch)
    TRIM(fastq_ch)
    NANOCHECK2(TRIM.out.trimmed_fastq)
}
process NANOCHECK1 {
    module 'nanoplot'
    debug true
    publishDir("results/nanoplot_before_trim", mode: 'copy')
    input: path sample_id
    output: path 'NanoPlot_CHECK_1', emit: nanoplot_files
    script:
    """
    NanoPlot -t 15 --fastq $sample_id --tsv_stats -o NanoPlot_CHECK_1
    """
}
process TRIM {
    module 'nanofilt'
    debug true
    publishDir("results/trimmed_fastq", mode: 'copy')
    input: path sample_id
    output: path 'sample_id.trimmed.fastq', emit: trimmed_fastq
    script:
    """
    NanoFilt -q 10 $sample_id > sample_id.trimmed.fastq
    """
}
process NANOCHECK2 {
    module 'nanoplot'
    debug true
    publishDir("results/nanoplot_after_trim", mode: 'copy')
    input: path sample_id
    output: path 'NanoPlot_CHECK_2', emit: nanoplot_files2
    script:
    """
    NanoPlot -t 15 --fastq $sample_id --tsv_stats -o NanoPlot_CHECK_2
    """
}

