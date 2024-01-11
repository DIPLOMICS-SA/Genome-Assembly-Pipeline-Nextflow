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
params.podsF = '/home/staukobong/pod5/'

pods_ch = Channel.fromPath(params.podsF, checkIfExists: true)
lineage_ch = Channel.fromPath(params.lineage, checkIfExists: true)

/*
 * Basecalling PORE5 files using Dorado
 */


process BASECALL {
    module 'dorado'
    debug true

    publishDir("${params.outdir}/bam_bascall", mode: 'copy')

    input:
    path sample_id

    output:
    path 'sample_id.bam' , emit: bamfiles_complete

    script:
    """
    dorado basecaller dna_r10.4.1_e8.2_400bps_hac@v4.2.0 $sample_id > sample_id.bam
    """
}


/*
 * Convert fastq files to bam files and concatenate the files
 */


process CONVERT {
    module 'samtools/1.9'
    debug true

    publishDir("${params.outdir}/coverted_fastq", mode: 'copy')

    input:
    path sample_id

    output:
    path 'sample_id.fastq', emit: fastq_files

    script:
    """
    samtools bam2fq $sample_id > sample_id.fastq
    """
}


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
    flye --nano-raw $sample_id -o assembly -i 1 -t 6
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
    racon -m 8 -x -8 -g -6 -w 500 -t 6 sample_id.trimmed.fastq sample_id.sam assembly/assembly.fasta > Racon_polished.fasta
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
    medaka_consensus -t 6 -m r1041_e82_400bps_hac_v4.2.0 -i sample_id.trimmed.fastq -d Racon_polished.fasta -o medaka_polished
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
    busco -m genome -in medaka_polished/consensus.fasta -o Busco_outputs -l eukaryota_odb10.2020-09-10.tar --metaeuk_parameters METAEUK_PARAMETERS --offline
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
    quast -t 6 -o Quast_output -r assembly/assembly.fasta Racon_polished.fasta medaka_polished/consensus.fasta
    """

}

 

/*
========================================================================================
                                Create default workflow
========================================================================================
*/

workflow {
    BASECALL(pods_ch)
    CONVERT(BASECALL.out.bamfiles_complete)
    FASTQC1(CONVERT.out.fastq_files)
    TRIM(CONVERT.out.fastq_files)
    FASTQC2(TRIM.out.trimmed_fastq)
    ASSEMBLY(TRIM.out.trimmed_fastq)
    MAPPINGS(TRIM.out.trimmed_fastq.combine(ASSEMBLY.out.Assembly_files))
    POLISH1(TRIM.out.trimmed_fastq.combine(MAPPINGS.out.Mapped_files.combine(ASSEMBLY.out.Assembly_files)))
    POLISHMED(TRIM.out.trimmed_fastq.combine(POLISH1.out.Polished_files))
    BUSCOstat(POLISH1.out.Polished_files.combine(lineage_ch))
    assemblyStats(ASSEMBLY.out.Assembly_files.combine(POLISH1.out.Polished_files.combine(POLISHMED.out.Polished_files2)))

}
