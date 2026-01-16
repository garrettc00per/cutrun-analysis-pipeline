#!/usr/bin/env nextflow

/*
 * BAM Processing Module
 * Quality filtering and fragment extraction for CUT&RUN data
 */

// Parameters specific to this module
params.outdir = "${projectDir}/results"
params.blacklist = "${projectDir}/ENCFF356LFX.bed.gz"
params.genome_sizes = "${projectDir}/GRCh38.p13.chrom.sizes"

process BAM_PROCESSING {
    publishDir "${params.outdir}/bam_processing", mode: 'copy'
    
    cpus 8
    
    input:
    tuple val(sample_id), path(bam_file)
    
    output:
    tuple val(sample_id), path ("${sample_id}_final.bam"), emit: bam
    tuple val(sample_id), path ("${sample_id}_final.bam.bai"), emit: bai
    tuple val(sample_id), path ("${sample_id}_final.clean.bedpe"), emit: bedpe
    tuple val(sample_id), path ("${sample_id}_final.fragments.bed"), emit: fragments
    tuple val(sample_id), path ("${sample_id}_final.fragments.sorted.bedgraph"), emit: bedgraph
    tuple val(sample_id), path ("${sample_id}_run.log"), emit: log
    
    script:
    """
    exec > ${sample_id}_run.log 2>&1
    set -euo pipefail
    
    echo "Processing ${sample_id}..."
    
    # Sort BAM by query name
    samtools sort -n -@ ${task.cpus} ${bam_file} -o ${sample_id}_sort.bam
    
    # Filter for high-quality, paired reads
    samtools view -F 3852 -f 2 -q 30 -h -@ ${task.cpus} -o ${sample_id}_sort_filt.bam ${sample_id}_sort.bam
    rm ${sample_id}_sort.bam
    
    # Remove blacklisted and non-chromosomal reads
    bedtools intersect -v -a ${sample_id}_sort_filt.bam -b ${params.blacklist} \\
      | samtools view -h -@ ${task.cpus} \\
      | awk -F '\\t' '!(\$3 ~ /chrM|chrUn|random|chrEBV/)' \\
      | samtools view -bh -@ ${task.cpus} > ${sample_id}_sort_filt2.bam
    rm ${sample_id}_sort_filt.bam
    
    # Name sort BAM file
    samtools sort -n -@ ${task.cpus} ${sample_id}_sort_filt2.bam -o ${sample_id}_final_namesort.bam
    rm ${sample_id}_sort_filt2.bam
    
    # Coordinate sort the same bam file
    samtools sort -@ ${task.cpus} ${sample_id}_final_namesort.bam -o ${sample_id}_final.bam
    
    # Index coordinate sorted bam file
    samtools index -@ ${task.cpus} ${sample_id}_final.bam
    
    # Convert to paired-end BEDPE format
    bedtools bamtobed -bedpe -i ${sample_id}_final_namesort.bam > ${sample_id}_final.bedpe
    rm ${sample_id}_final_namesort.bam
    
    # Filter pairs
    awk '\$1==\$4 && \$6-\$2 < 1000 {print \$0}' ${sample_id}_final.bedpe > ${sample_id}_final.clean.bedpe
    rm ${sample_id}_final.bedpe
    
    # Extract fragment coordinates
    cut -f 1,2,6 ${sample_id}_final.clean.bedpe | sort -k1,1 -k2,2n -k3,3n > ${sample_id}_final.fragments.bed
    
    # Generate bedGraph
    bedtools genomecov -bg -i ${sample_id}_final.fragments.bed -g ${params.genome_sizes} > ${sample_id}_final.fragments.bedgraph
    
    # Sort bedGraph
    sort -k1,1 -k2,2n ${sample_id}_final.fragments.bedgraph > ${sample_id}_final.fragments.sorted.bedgraph
    rm ${sample_id}_final.fragments.bedgraph
    
    echo "${sample_id} done."
    """
}

// Workflow that can be imported by main.nf
workflow RUN_BAM_PROCESSING {
    take:
    sample_ch
    
    main:
    BAM_PROCESSING(sample_ch)
    
    emit:
    bam = BAM_PROCESSING.out.bam
    bai = BAM_PROCESSING.out.bai
    bedgraph = BAM_PROCESSING.out.bedgraph
}
