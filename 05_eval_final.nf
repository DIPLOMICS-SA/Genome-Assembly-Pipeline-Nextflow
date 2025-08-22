nextflow.enable.dsl = 2

params.threads = params.threads ?: 15
params.lineage = params.lineage ?: 'eukaryota_odb10'
params.species_name = params.species_name ?: 'unknown_species'

process BUSCOstat_final {
    container '/home/apps/chpc/bio/busco/5.8.0/busco_v5.8.0_cv1.sif'
    publishDir "results/Busco_results", mode: 'copy'
    input:
    path assembly_file
    output:
    path 'Busco_output', emit: busco_results
    script:
    """
    busco -i ${assembly_file} \\
          -o Busco_output \\
          -l ${params.lineage} \\
          --offline \\
          -m genome \\
          --download_path ${workflow.projectDir}/busco_downloads \\
          -c ${params.threads}
    """
}

process QUAST_final {
    module 'quast/5.2.0'
    publishDir "results/quast_report", mode: 'copy'
    input:
    path assembly_file
    output:
    path 'Quast_result', emit: quast_results
    script:
    """
    quast.py -t ${params.threads} -o Quast_result --gene-finding --eukaryote ${assembly_file} --fragmented
    """
}

workflow {
    assembly_ch = Channel.fromPath(params.assembly)

    BUSCOstat_final(assembly_ch)
    QUAST_final(assembly_ch)
}
