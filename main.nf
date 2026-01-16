// main.nf

nextflow.enable.dsl=2

// Import modules
include { PEAK_CALLING } from './modules/peak_calling'

// Default params (can be overridden by config or command line)
params.bedgraph_dir = "${projectDir}/results/normalized_bedgraphs"
params.outdir = 'results'
params.seacr_script = "${projectDir}/SEACR_1.3.sh"
params.test_mode = false

workflow {
    // Load all normalized bedgraphs
    Channel
        .fromPath("${params.bedgraph_dir}/*_normalized.bedgraph")
        .map { file ->
            def filename = file.name
            // Parse: WT_SMARCB1_R1_normalized.bedgraph
            def parts = filename.replace('_normalized.bedgraph', '').tokenize('_')
            def genotype = parts[0]
            def target = parts[1]
            def replicate = parts[2]
            def sample_id = "${genotype}_${target}_${replicate}"
            [sample_id, genotype, target, replicate, file]
        }
        .filter { params.test_mode ? it[1] == 'WT' : true }  // WT only in test mode
        .set { bedgraph_ch }
    
    // Debug: see what we're loading
    bedgraph_ch.view { "Loading: ${it[0]}" }
    
    // Run peak calling
    PEAK_CALLING(bedgraph_ch)
    
    // View summaries
    PEAK_CALLING.out.filter_summary.view()
    PEAK_CALLING.out.overlap_summary.view()
}
