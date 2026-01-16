#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Import modules
include { RUN_BAM_PROCESSING as BAM_PROCESSING } from './modules/bam_processing'
include { BEDGRAPH_NORMALIZATION } from './modules/bedgraph_normalization'

// Parameters
params.sample_list = 'samples.txt'
params.bam_dir = '/home/ec2-user/cutnrun/full_run/bams'
params.norm_factors = 'normalization_factors.tsv'
params.outdir = 'results'

workflow {
    // Read sample list (just sample IDs, one per line)
    sample_ch = Channel
        .fromPath(params.sample_list)
        .splitText()
        .map { it.trim() }
        .map { sample_id -> [sample_id, file("${params.bam_dir}/${sample_id}.bam")] }
    
    // Run BAM processing
    BAM_PROCESSING(sample_ch)
    
    // Read normalization factors
    norm_factors_ch = Channel
        .fromPath(params.norm_factors)
        .splitCsv(sep: '\t', header: true)
        .map { row -> [row.Prefix, row.Factor.toFloat()] }
    
    // Run bedgraph normalization
    BEDGRAPH_NORMALIZATION(
        BAM_PROCESSING.out.bedgraph,
        norm_factors_ch
    )
}
