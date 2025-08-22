nextflow.enable.dsl = 2

params.threads = params.threads ?: 15
params.species_name = params.species_name ?: 'unknown_species'

process POLISH1 {
    module 'racon/1.5.0'
    publishDir "results/Racon_results", mode: 'copy'
    input:
    path fastq_file
    path sam_file
    path assembly_file
    output:
    path "${params.species_name}_Racon_polished.fasta", emit: polished_assembly
    script:
    """
    racon -m 3 -x -5 -g -4 -w 500 -t ${params.threads} ${fastq_file} ${sam_file} ${assembly_file} > ${params.species_name}_Racon_polished.fasta
    """
}

workflow {
    fastq_ch = Channel.fromPath(params.fastq)
    sam_ch = Channel.fromPath(params.sam)
    assembly_ch = Channel.fromPath(params.assembly)

    POLISH1(fastq_ch, sam_ch, assembly_ch)
}

