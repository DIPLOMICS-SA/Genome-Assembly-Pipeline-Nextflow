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
params.fastfiles = 'reads.fastq'

fastfiles_ch = Channel.fromPath(params.fastfiles, checkIfExists: true)


/*
 * Check quality of sequencing reads using FASTQC
 */

process FASTQC1 {
    debug true

    publishDir("${params.outdir}/fastqc_before_trim", mode: 'copy')

    input:
    path sample_id

    output:
    path 'NanoPlot_FASTQC_1', emit: fastqc_files

    script:
    """
    NanoPlot -t 15 --fastq $sample_id --tsv_stats -o NanoPlot_FASTQC_1
    """

}


/*
 * Trim fastq files after base calling using Nanofilt
 */

process TRIM {
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
 * Check quality of sequencing reads using FASTQC
 */

process FASTQC2 {
    debug true

    publishDir("${params.outdir}/fastqc_trim", mode: 'copy')

    input:
    path sample_id

    output:
    path 'NanoPlot_FASTQC_2', emit: fastqc_files2

    script:
    """
    NanoPlot -t 15 --fastq $sample_id --tsv_stats -o NanoPlot_FASTQC_2
    """

}


/*
 * Assemble the reads using FLYE
 */

process ASSEMBLY {
    debug true

    publishDir("${params.outdir}/Flye_results", mode: 'copy')

    input:
    path sample_id

    output:
    path 'assembly', emit: Assembly_files

    script:
    """
    flye --nano-raw $sample_id -o assembly --asm-coverage 40 -g 2.87g -t 15
    """

}


/*
 * Genome assembly assessment using Busco
 */

process BUSCOstat1 {
    debug true

    publishDir("${params.outdir}/Busco_results", mode: 'copy')

    input:
    path lineage

    output:
    path 'Busco_outputs1'

    script:
    """
    busco -i assembly/assembly.fasta -o Busco_outputs1 -m genome -l eukaryota_odb10 --metaeuk_parameters METAEUK_PARAMETERS --offline
    """

}


/*
 * Assembly evaluation using QUAST
 */

process assemblyStats1 {
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
    debug true

    publishDir("${params.outdir}/Racon_results", mode: 'copy')

    input:
    path sample_id

    output:
    path 'Racon_polished.fasta', emit: Polished_files

    script:
    """
    racon -m 8 -x -8 -g -6 -t 15 $sample_id sample_id.sam assembly/assembly.fasta > Racon_polished.fasta
    """

}


/*
 * Genome assembly assessment using Busco
 */

process BUSCOstat2 {
    debug true

    publishDir("${params.outdir}/Busco_results", mode: 'copy')

    input:
    path lineage

    output:
    path 'Busco_outputs2'

    script:
    """
    busco -i Racon_polished.fasta -o Busco_outputs2 -m genome -l eukaryota_odb10 --metaeuk_parameters METAEUK_PARAMETERS --offline
    """

}


/*
 * Assembly evaluation using QUAST
 */

process assemblyStats2 {
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

    FASTQC1(fastfiles_ch)
    TRIM(fastfiles_ch)
    FASTQC2(TRIM.out.trimmed_fastq)
    ASSEMBLY(TRIM.out.trimmed_fastq)
    BUSCOstat1(ASSEMBLY.out.Assembly_files)
    assemblyStats1(ASSEMBLY.out.Assembly_files)
    MAPPINGS(TRIM.out.trimmed_fastq.combine(ASSEMBLY.out.Assembly_files))
    POLISH1(fastfiles_ch.combine(MAPPINGS.out.Mapped_files.combine(ASSEMBLY.out.Assembly_files)))
    BUSCOstat2(POLISH1.out.Polished_files)
    assemblyStats2(POLISH1.out.Polished_files)

}
