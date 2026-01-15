// Parameters
params.outdir = 'results'

// modules/bedgraph_normalization.nf

process NORMALIZE_BEDGRAPH {
    tag "$sample_id"
    publishDir "${params.outdir}/normalized_bedgraphs", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bedgraph), val(norm_factor)
    
    output:
    tuple val(sample_id), path("${sample_id}_normalized.bedgraph")
    
    script:
    """
    awk -v factor="${norm_factor}" 'BEGIN{OFS="\\t"} {
        \$4 = \$4 * factor
        print
    }' ${bedgraph} > ${sample_id}_normalized.bedgraph
    """
}

process BASELINE_SUBTRACT {
    tag "$sample_id"
    publishDir "${params.outdir}/baseline_normalized", mode: 'copy'
    
    input:
    tuple val(sample_id), path(normalized_bg)
    
    output:
    tuple val(sample_id), path("${sample_id}_baseline.bedgraph")
    
    script:
    """
    # Find minimum value
    min=\$(awk 'NR==1 {min=\$4} \$4<min {min=\$4} END {print min}' ${normalized_bg})
    
    # Subtract minimum
    awk -v min=\$min 'BEGIN{OFS="\\t"} {
        corrected = \$4 - min
        if (corrected < 0) corrected = 0
        print \$1, \$2, \$3, corrected
    }' ${normalized_bg} > ${sample_id}_baseline.bedgraph
    """
}

process AVERAGE_REPLICATES {
    tag "${genotype}_${antibody}"
    publishDir "${params.outdir}/averaged_bedgraphs", mode: 'copy'
    
    input:
    tuple val(genotype), val(antibody), path(r1_bg), path(r2_bg)
    
    output:
    tuple val(genotype), val(antibody), path("${genotype}_${antibody}_averaged.bedgraph")
    
    script:
    """
    bedtools unionbedg -i ${r1_bg} ${r2_bg} | \\
    awk 'BEGIN{OFS="\\t"} {
        sum=0; count=0
        for(i=4; i<=NF; i++) {
            if(\$i!=".") {sum+=\$i; count++}
        }
        if(count>0) print \$1, \$2, \$3, sum/count
        else print \$1, \$2, \$3, "0"
    }' > ${genotype}_${antibody}_averaged.bedgraph
    """
}

workflow BEDGRAPH_NORMALIZATION {
    take:
    bedgraph_ch          // channel: [sample_id, bedgraph_file]
    norm_factors_ch      // channel: [sample_id, norm_factor]
    
    main:
    // Join bedgraphs with normalization factors
    bedgraph_with_factors = bedgraph_ch.join(norm_factors_ch)
    
    // Normalize
    normalized = NORMALIZE_BEDGRAPH(bedgraph_with_factors)
    
    // Baseline subtract
    baseline = BASELINE_SUBTRACT(normalized)
    
    // Parse sample IDs to extract genotype and antibody
    // Assumes format: GENOTYPE_ANTIBODY_REP
    baseline_parsed = baseline.map { sample_id, bedgraph ->
        def parts = sample_id.tokenize('_')
        def genotype = parts[0]
        def antibody = parts[1]
        def rep = parts[2]
        [genotype, antibody, rep, bedgraph]
    }
    
    // Group by genotype and antibody
    rep_pairs = baseline_parsed
        .groupTuple(by: [0, 1])
        .filter { it[3].size() == 2 }  // Only pairs with both replicates
        .map { genotype, antibody, reps, bedgraphs ->
            [genotype, antibody, bedgraphs[0], bedgraphs[1]]
        }
    
    // Average replicates
    averaged = AVERAGE_REPLICATES(rep_pairs)
    
    emit:
    normalized = normalized
    baseline = baseline
    averaged = averaged
}
