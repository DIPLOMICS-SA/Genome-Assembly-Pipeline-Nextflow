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


params.pipelineVersion = "1.0.0"

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
params.fastfiles = '/home/barcode12_fastq_pass.fastq.gz'

fastfiles_ch = Channel.fromPath(params.fastfiles, checkIfExists: true)

/*
 * Check quality of sequencing reads using FASTQC
 */

process FASTQC1 {
    module 'nanoplot'
    debug true

    publishDir("${params.outdir}/fastqc_before_trim", mode: 'copy')

    input:
    path sample_id

    output:
    path 'NanoPlot_FASTQC_1', emit: fastqc_files

    script:
    """
    NanoPlot -t 8 --fastq $sample_id --tsv_stats --plots hex dot -o NanoPlot_FASTQC_1
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
    NanoFilt -l 200 -q 15 --headcrop 50 --tailcrop 50  $sample_id > sample_id.trimmed.fastq
    """
}


/*
 * Check quality of sequencing reads using FASTQC
 */

process FASTQC2 {
    module 'nanoplot'
    debug true

    publishDir("${params.outdir}/fastqc_trim", mode: 'copy')

    input:
    path sample_id

    output:
    path 'NanoPlot_FASTQC_2', emit: fastqc_files2

    script:
    """
    NanoPlot -t 8 --fastq $sample_id --tsv_stats --plots hex dot -o NanoPlot_FASTQC_2
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
    flye --nano-raw $sample_id -o assembly --asm-coverage 40 -g 2.87g -t 8
    """

}


/*
 * Mapping the reads using minimap2
 */

process MAPPINGS {
    module 'minimap2'
    debug true

    publishDir("${params.outdir}/sam_file", mode: 'copy')

    input:
    path sample_id

    output:
    path 'sample_id.sam', emit: Mapped_files

    script:
    """
    minimap2 -ax map-ont -t 6 assembly/assembly.fasta sample_id.trimmed.fastq > sample_id.sam
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
    racon -m 8 -x -8 -g -6 -w 500 -t 4 sample_id.trimmed.fastq sample_id.sam assembly/assembly.fasta > Racon_polished.fasta
    """

}


/*
 * Polishing assembly using Medaka
 */

process POLISHMED {
    module 'medaka/1.11.1'
    debug true

    publishDir("${params.outdir}/Medaka_results", mode: 'copy')

    input:
    path sample_id

    output:
    path 'medaka_polished', emit: Polished_files2

    script:
    """
    medaka_consensus -t 4 -m r1041_e82_400bps_hac_v4.2.0 -i sample_id.trimmed.fastq -d Racon_polished.fasta -o medaka_polished
    """
    
}


/*
 * Genome assembly assessment using Busco
 */

process BUSCOstat {
    module 'busco'
    debug true

    publishDir("${params.outdir}/Busco_results", mode: 'copy')

    input:
    path lineage

    output:
    path 'Busco_outputs'

    script:
    """
    busco -i medaka_polished/consensus.fasta -o Busco_outputs -m genome -l eukaryota_odb10 --metaeuk_parameters METAEUK_PARAMETERS --offline
    """

}

/*
 * Assembly evaluation using QUAST
 */

process assemblyStats {
    module 'quast/4.6.3'
    debug true

    publishDir("${params.outdir}/quast_report", mode: 'copy')

    input:
    path sample_id

    output:
    path 'Quast_output'

    script:
    """
    quast -t 4 -o Quast_output -r assembly/assembly.fasta Racon_polished.fasta consensus.fasta
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
    MAPPINGS(TRIM.out.trimmed_fastq.combine(ASSEMBLY.out.Assembly_files))
    POLISH1(TRIM.out.trimmed_fastq.combine(MAPPINGS.out.Mapped_files.combine(ASSEMBLY.out.Assembly_files)))
    POLISHMED(TRIM.out.trimmed_fastq.combine(POLISH1.out.Polished_files))
    BUSCOstat(POLISHMED.out.Polished_files2.combine(lineage_ch))
    assemblyStats(ASSEMBLY.out.Assembly_files.combine(POLISH1.out.Polished_files.combine(POLISHMED.out.Polished_files2)))

}
