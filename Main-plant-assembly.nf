#!/usr/bin/env nextflow

/*
========================================================================================
                         plant_genome_ONT_assembly
========================================================================================

[Veldkos | Field Food | Spekboom Projects]
Contact: Larysha Rothmann <larysha@cengen.co.za>

*/

nextflow.enable.dsl = 2

if (params.help) {

    log.info """
    =============================================================================
                         plant_genome_ONT_assembly ~ version 1.0.0
    =============================================================================
   
   The downstream analysis and assembly of ONT basecalled reads, including read QC 
   (chopper and nanoplot), assembly (flye), assembly polishing (racon and medaka),
   assembly QC (quast, busco and blobtools), assembly assessment (jellyfish and genomescope)
   and removing contaminants (custom script). 
   
   NOTE: this script depends on `env.sh` for environment modules and `nextflow.config` for PBS settings.

    USAGE: 
    nextflow run ont_assembly.nf OPTIONS

    OPTIONS:
    <REQUIRED> -
    --reads : full path to the concatenated file of input ONT basecalled reads in fastq.gz format
    --outdir : full path to an output directory for results (default is ./results)

    <OPTIONAL> -
    --help : disply this message
    --busco_model : full path to busco model for core genes - default is viridiplantae_odb10 in 
                    "/mnt/lustre/users/lrothmann/busco_modules"
    --approx_size : the approximate or estimated genome size in MB (jellyfish) - default is 80000000                    

    NEXTFLOW OPTIONS:
    -resume : continue pipeline from where it left off after error / cluster time out
    -with-trace : tabular file with tracings of each processes [FILE] 

"""
    exit 0
}


/*
================================================================
            DEFINE PARAMETERS
================================================================
*/

params.version = "1.0.0"
params.reads = '/mnt/lustre/groups/CBBI1622/RESULTS/Polyphylla-test/oxalis_polyphylla_barcode02_fastq_pass.fastq'
params.outdir = "./" 
params.busco_model = "/mnt/lustre/groups/CBBI1622/LINEAGES/lineage_files"
params.approx_size = "80000000" 

/*
================================================================
            DEFINE PROCESSES
================================================================
*/

process nanofilt {
    // the input file/s are defined by the user as --reads in fastq.gz format
    // only reads with Q10 are retained and a read head/tail of 5bps is trimmed
    module 'nanofilt'
    debug true

    publishDir("${params.outdir}/nanofilt/", mode: 'copy')

    input:
    path raw_reads
    output:
    path 'pri_filt.fastq', emit: filtered_reads

    script:
    """
    NanoFilt -l 200 -q 15 --headcrop 10 --tailcrop 10 $raw_reads > pri_filt.fastq
    """
}

process nanoplot {
    // Using the output from the chopper process as input
    // assess the quality of the reads and basic stats
    module 'nanoplot'
    debug true

    publishDir("${params.outdir}/nanoplot", mode: 'copy')

    input:
    path filtered_reads
    output:
    path nanoplot, emit: nanoplot_out

    script:
    """
    NanoPlot \
        -t 8 \
        --color blue \
        --fastq $filtered_reads \
        --tsv_stats \
        -o nanoplot \
        --no_static
    """
}

process flye {
    // Using the output from the chopper process as input
    // output a draft primary genome assembly (collapsed)
    module 'flye/2.9'
    debug true

    publishDir("${params.outdir}/flye", mode: 'copy')

    input:
    path filtered_reads
    output:
    path 'flye', emit: assembly_all
    path "assembly.fasta", emit: assembly_out

    script:
    """
    flye \
        --nano-raw \
        $filtered_reads \
        --asm-coverage 40 \
        --keep-haplotypes \
        -g 0.66g \
        --threads 20 \
        --out-dir .
    """
}

process minimap {
    //  map the reads back to the assembly
    module 'minimap2/2.24:samtools/1.9 '
    debug true

    publishDir("${params.outdir}/minimap", mode: 'copy')

    input:
    path assembly_out
    path filtered_reads
    output:
    path 'aln_sorted.bam', emit: minimapped
    path 'aln_sorted.bam.bai', emit: minimapped2

    script:
    """
    minimap2 \
    -ax \
    map-ont $assembly_out \
    $filtered_reads \
    > aln.sam

    samtools view \
    -bS aln.sam \
    > aln.bam

    samtools sort \
    aln.bam \
    > aln_sorted.bam

    samtools index \
    aln_sorted.bam 
    """
}

process racon {
    //  polish the assembly with racon
    module 'racon'
    debug true

    publishDir("${params.outdir}/racon", mode: 'copy')

    input:
    path assembly_out
    path filtered_reads
    path minimapped
    output:
    path 'racon.fasta', emit: racon_out

    script:
    """
    samtools view -h -o aln_sorted.sam $minimapped

    racon \
    -q 15 \
    -m 8 \
    -x -6 \
    -g -8 \
    $filtered_reads \
    aln_sorted.sam \
    $assembly_out \
    > racon.fasta
    """
}

process quast {
    //  based on the output of medaka, run assembly QC
    module 'quast/5.2.0'
    debug true

    publishDir("${params.outdir}/quast", mode: 'copy')

    input:
    path racon_out
    output:
    path quast, emit: quast_qc

    script:
    """
    quast.py \
        -o quast \
        --gene-finding \
        --eukaryote \
        $racon_out
    """
}

process busco {
    //  based on the output of medaka, run assembly QC
    module 'busco/5.6.1'
    debug true

    publishDir("${params.outdir}/busco", mode: 'copy')

    input:
    path racon_out
    output:
    path 'busco', emit: busco_qc

    script:
    """
    busco \
        -i $racon_out \
        -l eukaryota_odb10 \
        -o busco \
        -f \
        -m genome \
        --offline \
        --download_path ${params.busco_model}
    """
}

process blast {
    //  blast the assembly to the NCBI nucleotide database 
    module 'blast'
    debug true

    publishDir("${params.outdir}/blast", mode: 'copy')

    input:
    path racon_out
    output:
    path 'blast_nt.out', emit: blast_nt

    script:
    """
    which blastn
    env | grep DB

    blastn \
        -db nt \
        -query $racon_out \
        -out blast_nt.out \
        -max_target_seqs 20 \
        -max_hsps 1 \
        -evalue 1e-35 \
        -outfmt "6 qseqid staxids bitscore std"
    """
}

process blobtools {
    //  blast the assembly to the NCBI nucleotide database
    module 'blobtools/1.1.1'
    debug true

    publishDir("${params.outdir}/blobout", mode: 'copy')

    input:
    path racon_out
    path minimapped
    path blast_nt
    output:
    path blobout, emit: blob

    script:
    """
    blobnames="/apps/chpc/bio/blobtools/1.1.1/data/names.dmp"
    blobnodes="/apps/chpc/bio/blobtools/1.1.1/data/nodes.dmp"

    blobtools create \
        -i $racon_out \
        -b $minimapped \
        --names \$blobnames \
        --nodes \$blobnodes \
        -t $blast_nt \
        -o blob_db

    blobtools view \
        -i blob_db.blobDB.json \
        -o blobout

    blobtools plot \
        -i blobdb_db.blobDB.json \
        -o blobout
    """
}

process jellyfish {
    module 'jellyfish/2.2.10'
    debug true

    publishDir("${params.outdir}/jelly", mode: 'copy')

    input:
    path filtered_reads
    output:
    path 'reads-kmers.histo', emit: jelly

    script:
    """
    jellyfish count \
    -C -m 21 -t 10 \
    -s ${params.approx_size} \
    $filtered_reads \
    -o reads-kmers.jf

    jellyfish histo \
    -t 10 reads-kmers.jf \
    > reads-kmers.histo
    """
}


/*
================================================================
            DEFINE WORKFLOW
================================================================
*/

workflow {
    def raw_reads = Channel.fromPath(params.reads, checkIfExists: true)
    nanofilt(raw_reads)
    nanofilt.out.filtered_reads.collect().set { filtered_reads }
    nanoplot(filtered_reads)
    flye(filtered_reads)
    minimap(flye.out.assembly_out, filtered_reads)
    racon(flye.out.assembly_out, filtered_reads, minimap.out.minimapped)
    quast(racon.out.racon_out)
    busco(racon.out.racon_out)
    blast(racon.out.racon_out)
    blobtools(blast.out.blast_nt, racon.out.racon_out, minimap.out.minimapped)
    jellyfish(filtered_reads)
}
