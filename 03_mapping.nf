nextflow.enable.dsl = 2

params.threads = params.threads ?: 15
params.species_name = params.species_name ?: 'unknown_species'

process MAPPINGS {
    publishDir "results/sam_file", mode: 'copy'
    input:
    path fastq_file
    path assembly_file
    output:
    path "${params.species_name}.sam", emit: mapped_sam
    script:
    """
    minimap2 -ax map-ont -t ${params.threads} ${assembly_file} ${fastq_file} > ${params.species_name}.sam
    """
}

workflow {
    MAPPINGS(params.fastq, params.assembly)
}

