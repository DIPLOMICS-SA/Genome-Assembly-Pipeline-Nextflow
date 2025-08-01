#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.assembly = null
params.lineage  = 'eukaryota_odb10'
params.outdir   = 'results'

workflow {
    Channel
        .fromPath(params.assembly, checkIfExists: true)
        .map { it -> [it] }         // wrap in tuple for DSL2
        .set { assembly_ch }

    BUSCOstat1(assembly_ch)
    assemblyStats1(assembly_ch)
}

/*
 * BUSCO evaluation before polishing
 */
process BUSCOstat1 {
    debug true
    container '/home/apps/chpc/bio/busco/5.8.0/busco_v5.8.0_cv1.sif'

    publishDir("${params.outdir}/Busco_results", mode: 'copy')

    input:
    tuple path(assembly_file)

    output:
    path 'Busco_outputs1'

    script:
    """
    busco -i ${assembly_file} \\
          -o Busco_outputs1 \\
          -l ${params.lineage} \\
          --metaeuk_parameters METAEUK_PARAMETERS \\
          --offline \\
          -m genome \\
          --download_path ${workflow.projectDir}/busco_downloads \\
          -c 12
    """
}

/*
 * QUAST evaluation before polishing
 */
process assemblyStats1 {
    module 'quast/5.2.0'
    debug true

    publishDir("${params.outdir}/quast_report", mode: 'copy')

    input:
    tuple path(assembly_file)

    output:
    path 'Quast_output1'

    script:
    """
    quast.py -t 15 -o Quast_output1 --gene-finding --eukaryote ${assembly_file} --fragmented
    """
}

