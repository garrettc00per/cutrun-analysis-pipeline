// modules/peak_calling.nf

process calculateIgGStats {
    tag "${genotype}_${replicate}"
    
    input:
    tuple val(genotype), val(replicate), path(igg_bedgraph)
    
    output:
    tuple val(genotype), val(replicate), path("${genotype}_IgG_${replicate}_stats.txt")
    
    script:
    """
    # Calculate mean and SD from IgG bedgraph (column 4 = intensity)
    awk '{sum+=\$4; sumsq+=\$4*\$4; n++} 
         END {
             mean=sum/n; 
             sd=sqrt((sumsq-sum*sum/n)/n);
             snr_threshold=mean+(${params.snr_multiplier}*sd);
             print "mean\\t"mean;
             print "sd\\t"sd;
             print "snr_threshold\\t"snr_threshold
         }' ${igg_bedgraph} > ${genotype}_IgG_${replicate}_stats.txt
    
    echo "IgG stats for ${genotype} ${replicate}:"
    cat ${genotype}_IgG_${replicate}_stats.txt
    """
}

process callPeaks {
    tag "${sample_id}"
    publishDir "${params.outdir}/seacr_peaks", mode: 'copy'
    
    cpus 2
    
    input:
    tuple val(sample_id), val(genotype), val(target), val(replicate), 
          path(sample_bedgraph), path(control_bedgraph)
    
    output:
    tuple val(sample_id), val(genotype), val(target), val(replicate),
          path("${sample_id}_SEACR_peaks.stringent.bed")
    
    script:
    """
    # Run SEACR with IgG control
    bash ${params.seacr_script} \
        ${sample_bedgraph} \
        ${control_bedgraph} \
        non \
        stringent \
        ${sample_id}_SEACR_peaks
    """
}

process filterPeaks {
    tag "${sample_id}"
    publishDir "${params.outdir}/filtered_peaks", mode: 'copy'
    
    input:
    tuple val(sample_id), val(genotype), val(target), val(replicate),
          path(peak_bed), path(stats_file)
    
    output:
    tuple val(sample_id), val(genotype), val(target), val(replicate),
          path("${sample_id}_filtered.bed"), 
          path("${sample_id}_filter_summary.txt")
    
    script:
    """
    # Extract SNR threshold from stats file
    SNR_THRESHOLD=\$(grep "^snr_threshold" ${stats_file} | cut -f2)
    FIXED_THRESHOLD=${params.fixed_threshold}
    FILTER_METHOD="${params.filter_method}"
    
    # Determine final threshold based on method
    if [ "\$FILTER_METHOD" == "snr" ]; then
        THRESHOLD=\$SNR_THRESHOLD
        METHOD_DESC="SNR only (mean + ${params.snr_multiplier}×SD)"
    elif [ "\$FILTER_METHOD" == "fixed" ]; then
        THRESHOLD=\$FIXED_THRESHOLD
        METHOD_DESC="Fixed threshold"
    elif [ "\$FILTER_METHOD" == "hybrid" ]; then
        THRESHOLD=\$(awk -v snr="\$SNR_THRESHOLD" -v fixed="\$FIXED_THRESHOLD" 'BEGIN{print (snr > fixed) ? snr : fixed}')
        METHOD_DESC="Hybrid (max of SNR and fixed)"
    else
        echo "Error: Unknown filter_method '\$FILTER_METHOD'"
        exit 1
    fi
    
    echo "Sample: ${sample_id}" > ${sample_id}_filter_summary.txt
    echo "Filter method: \$METHOD_DESC" >> ${sample_id}_filter_summary.txt
    echo "SNR threshold: \$SNR_THRESHOLD" >> ${sample_id}_filter_summary.txt
    echo "Fixed threshold: \$FIXED_THRESHOLD" >> ${sample_id}_filter_summary.txt
    echo "Applied threshold: \$THRESHOLD" >> ${sample_id}_filter_summary.txt
    echo "Original peaks: \$(wc -l < ${peak_bed})" >> ${sample_id}_filter_summary.txt
    
    # Filter peaks by max intensity (column 5) > threshold
    awk -v thresh="\$THRESHOLD" '\$5 > thresh' ${peak_bed} > ${sample_id}_filtered.bed
    
    echo "Filtered peaks: \$(wc -l < ${sample_id}_filtered.bed)" >> ${sample_id}_filter_summary.txt
    echo "Percent retained: \$(awk -v orig=\$(wc -l < ${peak_bed}) -v filt=\$(wc -l < ${sample_id}_filtered.bed) 'BEGIN{printf "%.1f", (filt/orig)*100}')%" >> ${sample_id}_filter_summary.txt
    
    cat ${sample_id}_filter_summary.txt
    """
}

process summarizePeakCalling {
    publishDir "${params.outdir}/peak_summary", mode: 'copy'
    
    input:
    path(all_summaries)
    
    output:
    path("peak_calling_summary.txt")
    
    script:
    """
    echo -e "Sample\\tOriginal\\tFiltered\\tPercent\\tSNR_Threshold\\tFixed_Threshold\\tApplied_Threshold\\tMethod" > peak_calling_summary.txt
    
    for summary in ${all_summaries}; do
        sample=\$(grep "^Sample:" \$summary | cut -d' ' -f2)
        method=\$(grep "^Filter method:" \$summary | cut -d':' -f2- | xargs)
        snr_thresh=\$(grep "^SNR threshold:" \$summary | awk '{print \$3}')
        fixed_thresh=\$(grep "^Fixed threshold:" \$summary | awk '{print \$3}')
        applied_thresh=\$(grep "^Applied threshold:" \$summary | awk '{print \$3}')
        orig=\$(grep "^Original peaks:" \$summary | awk '{print \$3}')
        filt=\$(grep "^Filtered peaks:" \$summary | awk '{print \$3}')
        pct=\$(grep "^Percent retained:" \$summary | awk '{print \$3}')
        echo -e "\$sample\\t\$orig\\t\$filt\\t\$pct\\t\$snr_thresh\\t\$fixed_thresh\\t\$applied_thresh\\t\$method"
    done >> peak_calling_summary.txt
    
    echo ""
    echo "=== Peak Calling Summary ==="
    column -t peak_calling_summary.txt
    """
}

process findReplicateOverlaps {
    tag "${genotype}_${target}"
    publishDir "${params.outdir}/replicate_overlaps", mode: 'copy'
    
    input:
    tuple val(genotype), val(target), path(r1_bed), path(r2_bed)
    
    output:
    tuple val(genotype), val(target), 
          path("${genotype}_${target}_R1_vs_R2_overlap.bed"),
          path("${genotype}_${target}_overlap_stats.txt")
    
    script:
    """
    # Find overlapping peaks between replicates
    bedtools intersect -a ${r1_bed} -b ${r2_bed} -wa -wb > ${genotype}_${target}_R1_vs_R2_overlap.bed
    
    # Generate statistics
    R1_COUNT=\$(wc -l < ${r1_bed})
    R2_COUNT=\$(wc -l < ${r2_bed})
    OVERLAP_COUNT=\$(wc -l < ${genotype}_${target}_R1_vs_R2_overlap.bed)
    
    echo "Sample: ${genotype}_${target}" > ${genotype}_${target}_overlap_stats.txt
    echo "R1 peaks: \$R1_COUNT" >> ${genotype}_${target}_overlap_stats.txt
    echo "R2 peaks: \$R2_COUNT" >> ${genotype}_${target}_overlap_stats.txt
    echo "Overlapping peaks: \$OVERLAP_COUNT" >> ${genotype}_${target}_overlap_stats.txt
    
    # Calculate overlap percentages
    R1_PCT=\$(awk -v overlap=\$OVERLAP_COUNT -v r1=\$R1_COUNT 'BEGIN{printf "%.1f", (overlap/r1)*100}')
    R2_PCT=\$(awk -v overlap=\$OVERLAP_COUNT -v r2=\$R2_COUNT 'BEGIN{printf "%.1f", (overlap/r2)*100}')
    
    echo "R1 overlap %: \$R1_PCT" >> ${genotype}_${target}_overlap_stats.txt
    echo "R2 overlap %: \$R2_PCT" >> ${genotype}_${target}_overlap_stats.txt
    
    echo "✅ ${genotype}_${target}: \$OVERLAP_COUNT overlapping peaks (\$R1_PCT% of R1, \$R2_PCT% of R2)"
    """
}

process summarizeReplicateOverlaps {
    publishDir "${params.outdir}/replicate_overlaps", mode: 'copy'
    
    input:
    path(all_stats)
    
    output:
    path("replicate_overlap_summary.txt")
    
    script:
    """
    echo -e "Genotype\\tTarget\\tR1_Peaks\\tR2_Peaks\\tOverlapping\\tR1_Overlap_%\\tR2_Overlap_%" > replicate_overlap_summary.txt
    
    for stats in ${all_stats}; do
        sample=\$(grep "^Sample:" \$stats | cut -d' ' -f2)
        r1=\$(grep "^R1 peaks:" \$stats | awk '{print \$3}')
        r2=\$(grep "^R2 peaks:" \$stats | awk '{print \$3}')
        overlap=\$(grep "^Overlapping peaks:" \$stats | awk '{print \$3}')
        r1_pct=\$(grep "^R1 overlap %:" \$stats | awk '{print \$4}')
        r2_pct=\$(grep "^R2 overlap %:" \$stats | awk '{print \$4}')
        
        # Split sample into genotype and target
        genotype=\$(echo \$sample | cut -d'_' -f1)
        target=\$(echo \$sample | cut -d'_' -f2)
        
        echo -e "\$genotype\\t\$target\\t\$r1\\t\$r2\\t\$overlap\\t\$r1_pct\\t\$r2_pct"
    done >> replicate_overlap_summary.txt
    
    echo ""
    echo "=== Replicate Overlap Summary ==="
    column -t replicate_overlap_summary.txt
    """
}

workflow PEAK_CALLING {
    take:
    bedgraph_channel  // [sample_id, genotype, target, replicate, bedgraph_path]
    
    main:
    // Separate IgG controls from targets
    igg_controls = bedgraph_channel
        .filter { it[2] == 'IgG' }
        .map { [it[1], it[3], it[4]] }  // [genotype, replicate, bedgraph]
    
    target_samples = bedgraph_channel
        .filter { it[2] != 'IgG' }
    
    // Calculate IgG statistics
    igg_stats = calculateIgGStats(igg_controls)
    
    // Match each target sample with its IgG control
    samples_with_controls = target_samples
        .map { [it[1], it[3], it[0], it[2], it[4]] }  // [genotype, replicate, sample_id, target, bedgraph]
        .combine(igg_controls.map { [it[0], it[1], it[2]] }, by: [0,1])  // Match on genotype+replicate
        .map { [it[2], it[0], it[3], it[1], it[4], it[5]] }  // [sample_id, genotype, target, replicate, sample_bg, control_bg]
    
    // Call peaks
    peaks = callPeaks(samples_with_controls)
    
    // Match peaks with their IgG stats for filtering
    peaks_with_stats = peaks
        .map { [it[1], it[3], it[0], it[2], it[4]] }  // [genotype, replicate, sample_id, target, peak_bed]
        .combine(igg_stats.map { [it[0], it[1], it[2]] }, by: [0,1])  // Match on genotype+replicate
        .map { [it[2], it[0], it[3], it[1], it[4], it[5]] }  // [sample_id, genotype, target, replicate, peak_bed, stats]
    
    // Filter peaks by threshold
    filtered = filterPeaks(peaks_with_stats)
    
    // Summarize filtering results
    all_filter_summaries = filtered.map { it[5] }.collect()
    filter_summary = summarizePeakCalling(all_filter_summaries)
    
    // Pair replicates for overlap analysis
    // filtered: [sample_id, genotype, target, replicate, filtered_bed, summary]
    r1_peaks = filtered
        .filter { it[3] == 'R1' }
        .map { [it[1], it[2], it[4]] }  // [genotype, target, r1_bed]
    
    r2_peaks = filtered
        .filter { it[3] == 'R2' }
        .map { [it[1], it[2], it[4]] }  // [genotype, target, r2_bed]
    
    // Combine R1 and R2 by genotype+target
    replicate_pairs = r1_peaks
        .join(r2_peaks, by: [0,1])  // Join on genotype and target
        // Results in: [genotype, target, r1_bed, r2_bed]
    
    // Find overlaps between replicates
    overlaps = findReplicateOverlaps(replicate_pairs)
    
    // Summarize overlaps
    all_overlap_stats = overlaps.map { it[3] }.collect()
    overlap_summary = summarizeReplicateOverlaps(all_overlap_stats)
    
    emit:
    filtered_peaks = filtered.map { [it[0], it[1], it[2], it[3], it[4]] }
    overlapping_peaks = overlaps.map { [it[0], it[1], it[2]] }  // [genotype, target, overlap_bed]
    filter_summary = filter_summary
    overlap_summary = overlap_summary
}
