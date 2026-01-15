#!/usr/bin/env nextflow

/*
 * CUT&RUN Analysis Pipeline
 * Main workflow orchestrator
 */

// Parameters
params.sample_list = "${projectDir}/samples.txt"
params.bam_dir = "/home/ec2-user/cutnrun/full_run/bams"
params.outdir = "${projectDir}/results"

// Print header
log.info """\
    CUT&RUN ANALYSIS PIPELINE
    =========================
    sample list    : ${params.sample_list}
    BAM directory  : ${params.bam_dir}
    output dir     : ${params.outdir}
    """
    .stripIndent()

// Include modules
include { BAM_PROCESSING } from './modules/bam_processing'

// Main workflow
workflow {
    // Create sample channel
    Channel
        .fromPath(params.sample_list)
        .splitText()
        .map { it.trim() }
        .map { sample_id -> 
            def bam_file = file("${params.bam_dir}/${sample_id}.bam")
            if (!bam_file.exists()) {
                error "BAM file not found: ${bam_file}"
            }
            tuple(sample_id, bam_file)
        }
        .set { samples_ch }
    
    // Run BAM processing
    BAM_PROCESSING(samples_ch)
}
