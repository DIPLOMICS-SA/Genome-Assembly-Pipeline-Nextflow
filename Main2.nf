#!/usr/bin/env nextflow

nextflow.enable.dsl=2


/*
 *
========================================================================================
         GENASS: Genome Assembly Pipeline for Nanopore Sequencing Data 
========================================================================================
 
# Homepage / Documentation
 GitHub - DIPLOMICS [1KSA Genome Assembly project]
 # Authors
 Setshaba Taukobong <setshaba.taukobong@diplomics.org.za> <sc.taukobong@gmail.com>

---------------------------------------------------------------------------------------
 *
 */


params.pipelineVersion = "3.0.0"

def helpMessage(){
 log.info"""
 =========================================================
      GENOME ASSEMBLY ~ version ${params.pipelineVersion}
 =========================================================
   Usage:

   Command for running the pipeline is as follows:

   nextflow run GenomeAssemblyPipeline.nf OPTIONS

   OPTIONS:

    NextFlow options [OPTIONAL]:
        Produce an html report with useful metrics of the pipeline [FILE]          -with-report
        Produce a tabular file with tracings of each processes [FILE]              -with-trace
        Produce an html graphic of all process executed [FILE]                     -with-timeline
        Resume running pipeline where it has left off after an error for example   -resume


"""}


/*
========================================================================================
                        Define parameters, channels and processes
========================================================================================
*/


/*
 * Define the default parameters
 */ 

params.outdir = 'results'
params.fastfiles = 'species_name_fastq_pass_con.fastq'

fastfiles_ch = Channel.fromPath(params.fastfiles, checkIfExists: true)


/*
 * Check quality of sequencing reads using NANOPLOT
 */

process NANOCHECK1 {
    module 'nanoplot'
    debug true

    publishDir("${params.outdir}/nanoplot_before_trim", mode: 'copy')

    input:
    path sample_id

    output:
    path 'NanoPlot_CHECK_1', emit: nanoplot_files

    script:
    """
    NanoPlot -t 15 --fastq $sample_id --tsv_stats -o NanoPlot_CHECK_1
    """

}


/*
 * Trim fastq files after base calling using Nanofilt
 */

process TRIM {
    module 'nanofilt'
    debug true

    publishDir("${params.outdir}/trimmed_fastq", mode: 'copy')

    input:
    path sample_id

    output:
    path 'sample_id.trimmed.fastq', emit: trimmed_fastq

    script:
    """
    NanoFilt -q 10 $sample_id > sample_id.trimmed.fastq
    """
}


/*
 * Check quality of sequencing reads using NANOPLOT after trimming
 */

process NANOCHECK2 {
    module 'nanoplot'
    debug true

    publishDir("${params.outdir}/nanoplot_after_trim", mode: 'copy')

    input:
    path sample_id

    output:
    path 'NanoPlot_CHECK_2', emit: nanoplot_files2

    script:
    """
    NanoPlot -t 15 --fastq $sample_id --tsv_stats -o NanoPlot_CHECK_2
    """

}


/*
 * Assemble the reads using FLYE
 */

process ASSEMBLY {
    module 'flye/2.9'
    debug true

    publishDir("${params.outdir}/Flye_results", mode: 'copy')

    input:
    path sample_id

    output:
    path 'assembly', emit: Assembly_files

    script:
    """
    flye --nano-raw $sample_id -o assembly --asm-coverage x -g x g -t 15
    """

}


/*
 * Genome assembly assessment before polishing using Busco
 */

process BUSCOstat1 {
    module 'busco/5.4.5'
    debug true

    publishDir("${params.outdir}/Busco_results", mode: 'copy')

    input:
    path lineage

    output:
    path 'Busco_outputs1'

    script:
    """
    busco -m genome -i assembly/assembly.fasta -o Busco_outputs1 -l eukaryota_odb10 --metaeuk_parameters METAEUK_PARAMETERS --offline
    """

}


/*
 * Assembly evaluation before polishing using QUAST
 */

process assemblyStats1 {
    module 'quast/4.6.3'
    debug true

    publishDir("${params.outdir}/quast_report", mode: 'copy')

    input:
    path sample_id

    output:
    path 'Quast_output1'

    script:
    """
    quast.py -t 15 -o Quast_output1 --gene-finding --eukaryote assembly/assembly.fasta --fragmented
    """

}


/*
 * Mapping the reads using minimap2
 */

process MAPPINGS {
    debug true

    publishDir("${params.outdir}/sam_file", mode: 'copy')

    input:
    path sample_id

    output:
    path 'sample_id.sam', emit: Mapped_files

    script:
    """
    minimap2 -ax map-ont -t 15 assembly/assembly.fasta sample_id.trimmed.fastq > sample_id.sam
    """

}

/*
 * Polishing assembly using Racon
 */

process POLISH1 {
    module 'racon/1.5.0'
    debug true

    publishDir("${params.outdir}/Racon_results", mode: 'copy')

    input:
    path sample_id

    output:
    path 'Racon_polished.fasta', emit: Polished_files

    script:
    """
    racon -m 3 -x -5 -g -4 -w 500 -t 15 $sample_id sample_id.sam assembly/assembly.fasta > Racon_polished.fasta
    """

}

/*
 * Genome assembly assessment after polishing using Busco
 */

process BUSCOstat2 {
    module 'busco/5.4.5'
    debug true

    publishDir("${params.outdir}/Busco_results", mode: 'copy')

    input:
    path lineage

    output:
    path 'Busco_outputs2'

    script:
    """
    busco -m genome -i Racon_polished.fasta -o Busco_outputs2 -l eukaryota_odb10 --metaeuk_parameters METAEUK_PARAMETERS --offline
    """

}


/*
 * Assembly evaluation after polishing using QUAST
 */

process assemblyStats2 {
    module 'quast/4.6.3'
    debug true

    publishDir("${params.outdir}/quast_report", mode: 'copy')

    input:
    path sample_id

    output:
    path 'Quast_output2'

    script:
    """
    quast.py -t 15 -o Quast_output2 --gene-finding --eukaryote Racon_polished.fasta --fragmented
    """

}


 
/*
========================================================================================
                                Create default workflow
========================================================================================
*/

workflow {
    NANOCHECK1(fastfiles_ch)
    Trim(fastfiles_ch)
    NANOCHECK2(TRIM.out.trimmed_fastq)
    ASSEMBLY(TRIM.out.trimmed_fastq)
    BUSCOstat1(ASSEMBLY.out.Assembly_files)
    assemblyStats1(ASSEMBLY.out.Assembly_files)
    MAPPINGS(TRIM.out.trimmed_fastq.combine(ASSEMBLY.out.Assembly_files))
    POLISH1(fastfiles_ch.combine(MAPPINGS.out.Mapped_files.combine(ASSEMBLY.out.Assembly_files)))
    BUSCOstat2(POLISH1.out.Polished_files2)
    assemblyStats2(POLISH1.out.Polished_files2)

}
