// Parameters
params.outdir = 'results'

process CREATE_BIGWIG {
    tag "$sample_id"
    publishDir "${params.outdir}/bigwigs_scaled", mode: 'copy'
    
    cpus 4
    
    input:
    tuple val(sample_id), path(bam), path(bai), val(scale_factor)
    
    output:
    tuple val(sample_id), path("${sample_id}_scaled.bw")
    
    script:
    """
    bamCoverage \
        -b ${bam} \
        -o ${sample_id}_scaled.bw \
        --scaleFactor ${scale_factor} \
        --normalizeUsing None \
        --binSize 10 \
        --ignoreDuplicates \
        --extendReads \
        -p ${task.cpus}
    """
}

process SUBTRACT_IGG {
    tag "${sample_id}"
    publishDir "${params.outdir}/bigwigs_igg_subtracted", mode: 'copy'
    
    cpus 4
    
    input:
    tuple val(sample_id), path(chip_bw), path(igg_bw)
    
    output:
    tuple val(sample_id), path("${sample_id}_IgGsubtracted.bw")
    
    script:
    """
    bigwigCompare \
        -b1 ${chip_bw} \
        -b2 ${igg_bw} \
        --operation subtract \
        --binSize 50 \
        -p ${task.cpus} \
        --skipNAs \
        -o ${sample_id}_IgGsubtracted.bw
    """
}

process AVERAGE_REPLICATES {
    tag "${genotype}_${antibody}"
    publishDir "${params.outdir}/bigwigs_averaged", mode: 'copy'
    
    cpus 8
    
    input:
    tuple val(genotype), val(antibody), path(r1_bw), path(r2_bw)
    
    output:
    tuple val(genotype), val(antibody), path("${genotype}_${antibody}_avg50bp.bw")
    
    script:
    """
    bigwigCompare \
        -b1 ${r1_bw} \
        -b2 ${r2_bw} \
        -o ${genotype}_${antibody}_avg50bp.bw \
        --operation mean \
        --binSize 50 \
        -p ${task.cpus}
    """
}

workflow BIGWIG_GENERATION {
    take:
    bam_ch              // channel: [sample_id, bam, bai]
    norm_factors_ch     // channel: [sample_id, factor]
    
    main:
    // Join BAMs with normalization factors
    bam_with_factors = bam_ch.join(norm_factors_ch)
    
    // Create scaled bigwigs
    scaled_bw = CREATE_BIGWIG(bam_with_factors)
    
    // Parse sample IDs to extract components
    // Format: GENOTYPE_ANTIBODY_REP
    scaled_parsed = scaled_bw.map { sample_id, bw ->
        def parts = sample_id.tokenize('_')
        def genotype = parts[0]
        def antibody = parts[1]
        def rep = parts[2]
        [genotype, antibody, rep, sample_id, bw]
    }
    
    // Separate IgG controls and target samples
    igg_samples = scaled_parsed.filter { genotype, antibody, rep, sample_id, bw ->
        antibody == 'IgG'
    }
    
    target_samples = scaled_parsed.filter { genotype, antibody, rep, sample_id, bw ->
        antibody != 'IgG'
    }
    
    // Match each target with its IgG control using combine + filter
    // This allows one-to-many matching (one IgG -> multiple targets)
    targets_with_igg = target_samples
        .combine(igg_samples)
        .filter { target_genotype, target_antibody, target_rep, target_sample_id, target_bw,
                  igg_genotype, igg_antibody, igg_rep, igg_sample_id, igg_bw ->
            target_genotype == igg_genotype && target_rep == igg_rep
        }
        .map { target_genotype, target_antibody, target_rep, target_sample_id, target_bw,
               igg_genotype, igg_antibody, igg_rep, igg_sample_id, igg_bw ->
            [target_sample_id, target_bw, igg_bw]
        }
    
    // Subtract IgG
    igg_subtracted = SUBTRACT_IGG(targets_with_igg)
    
    // Parse for averaging
    igg_subtracted_parsed = igg_subtracted.map { sample_id, bw ->
        def parts = sample_id.tokenize('_')
        def genotype = parts[0]
        def antibody = parts[1]
        def rep = parts[2]
        [genotype, antibody, rep, bw]
    }
    
    // Group by genotype and antibody to pair replicates
    rep_pairs = igg_subtracted_parsed
        .groupTuple(by: [0, 1])
        .filter { it[3].size() == 2 }  // Only pairs with both replicates
        .map { genotype, antibody, reps, bws ->
            [genotype, antibody, bws[0], bws[1]]
        }
    
    // Average replicates
    averaged = AVERAGE_REPLICATES(rep_pairs)
    
    emit:
    scaled = scaled_bw
    igg_subtracted = igg_subtracted
    averaged = averaged
}
