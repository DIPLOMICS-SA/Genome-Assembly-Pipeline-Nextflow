nextflow.enable.dsl=2

params.fastfiles = null
params.assembly  = null
params.lineage   = 'eukaryota_odb10'
params.outdir    = 'results'

workflow {
    // Load input files
    Channel.fromPath(params.fastfiles, checkIfExists: true)
        .map { [it] }
        .set { reads_ch }

    Channel.fromPath(params.assembly, checkIfExists: true)
        .map { [it] }
        .set { asm_ch }

    // 1. Map reads to assembly
    MAPPINGS(reads_ch.combine(asm_ch))

    // 2. Polish using Racon
    POLISH1(
        reads_ch.combine(
            MAPPINGS.out.Mapped_files.combine(asm_ch)
        )
    )

    // 3. Evaluate polished assembly
    BUSCOstat2(POLISH1.out.Polished_files.map { [it] })
    assemblyStats2(POLISH1.out.Polished_files.map { [it] })
}

process MAPPINGS {
    debug true
    publishDir("${params.outdir}/sam_file", mode: 'copy')

    input:
    tuple path(fastq_file), path(assembly_file)

    output:
    path 'sample_id.sam', emit: Mapped_files

    script:
    """
    minimap2 -ax map-ont -t 15 ${assembly_file} ${fastq_file} > sample_id.sam
    """
}

process POLISH1 {
    module 'racon/1.5.0'
    debug true
    publishDir("${params.outdir}/Racon_results", mode: 'copy')

    input:
    tuple path(fastq_file), path(sam_file), path(assembly_file)

    output:
    path 'Racon_polished.fasta', emit: Polished_files

    script:
    """
    racon -m 3 -x -5 -g -4 -w 500 -t 15 ${fastq_file} ${sam_file} ${assembly_file} > Racon_polished.fasta
    """
}

process BUSCOstat2 {
    debug true
    container '/home/apps/chpc/bio/busco/5.8.0/busco_v5.8.0_cv1.sif'
    publishDir("${params.outdir}/Busco_results", mode: 'copy')

    input:
    tuple path(polished_assembly_file)

    output:
    path 'Busco_outputs2'

    script:
    """
    busco -i ${polished_assembly_file} \\
          -o Busco_outputs2 \\
          -l ${params.lineage} \\
          --metaeuk_parameters METAEUK_PARAMETERS \\
          --offline \\
          -m genome \\
          --download_path ${workflow.projectDir}/busco_downloads \\
          -c 12
    """
}

process assemblyStats2 {
    module 'quast/5.2.0'
    debug true
    publishDir("${params.outdir}/quast_report", mode: 'copy')

    input:
    tuple path(polished_assembly_file)

    output:
    path 'Quast_output2'

    script:
    """
    quast.py -t 16 -o Quast_output2 --gene-finding --eukaryote ${polished_assembly_file} --fragmented
    """
}

