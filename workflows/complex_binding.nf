// workflows/complex_binding.nf
// Generic workflow: Create consensus binding regions across multiple targets and genotypes
//
// USAGE EXAMPLES:
//
//   SWI/SNF complex (WT vs I315I):
//     nextflow run workflows/complex_binding.nf \
//       --genotypes WT,I315I \
//       --targets_to_merge SMARCA4,SMARCB1,SMARCE1
//
//   Compare mutants for single protein:
//     nextflow run workflows/complex_binding.nf \
//       --genotypes W281P,W281X \
//       --targets_to_merge SMARCB1
//
//   Histone mark co-occurrence:
//     nextflow run workflows/complex_binding.nf \
//       --genotypes WT,I315I \
//       --targets_to_merge H3K27me3,H3K4me3 \
//       --merge_distance 500
//
//   All parameters:
//     nextflow run workflows/complex_binding.nf \
//       --overlap_dir results/replicate_overlaps \
//       --output_dir results/my_analysis \
//       --genotypes WT,I315R \
//       --targets_to_merge SMARCA4,SMARCE1 \
//       --merge_distance 150

nextflow.enable.dsl=2

// Parameters
params.overlap_dir = 'results/replicate_overlaps'
params.output_dir = 'results/complex_binding'
params.genotypes = 'WT,I315I'
params.targets_to_merge = 'SMARCA4,SMARCB1,SMARCE1'
params.merge_distance = 100

// Parse comma-separated strings into lists
def genotypes = params.genotypes instanceof String ? params.genotypes.tokenize(',') : params.genotypes
def targets = params.targets_to_merge instanceof String ? params.targets_to_merge.tokenize(',') : params.targets_to_merge

// Validate parameters
if (genotypes.size() < 2) {
    error "ERROR: Must specify at least 2 genotypes for consensus analysis"
}
if (targets.size() < 1) {
    error "ERROR: Must specify at least 1 target to analyze"
}

log.info """
=========================================
Complex Binding Region Analysis
=========================================
Genotypes to compare: ${genotypes.join(', ')}
Targets to merge:     ${targets.join(', ')}
Merge distance:       ${params.merge_distance}bp
Input directory:      ${params.overlap_dir}
Output directory:     ${params.output_dir}
=========================================
"""

process createReproduciblePeaks {
    tag "${genotype}_${target}"
    publishDir "${params.output_dir}/reproducible_peaks", mode: 'copy'
    
    input:
    tuple val(genotype), val(target), path(overlap_bed)
    
    output:
    tuple val(genotype), val(target), path("${genotype}_${target}_reproducible_peaks.bed")
    
    script:
    """
    # Merge overlapping R1 and R2 peaks by taking union of coordinates
    awk '{
        chr = \$1;
        start = (\$2 < \$7 ? \$2 : \$7);
        end   = (\$3 > \$8 ? \$3 : \$8);
        print chr "\\t" start "\\t" end
    }' ${overlap_bed} > ${genotype}_${target}_reproducible_peaks.bed
    
    PEAK_COUNT=\$(wc -l < ${genotype}_${target}_reproducible_peaks.bed)
    echo "✅ ${genotype}_${target}: \$PEAK_COUNT reproducible peaks"
    """
}

process mergeTargets {
    tag "${genotype}"
    publishDir "${params.output_dir}/unified_regions", mode: 'copy'
    
    input:
    tuple val(genotype), path(peak_beds)
    
    output:
    tuple val(genotype), path("${genotype}_unified_regions.bed")
    
    script:
    """
    # Concatenate all target peaks
    cat ${peak_beds} > ${genotype}_all_peaks.bed
    
    # Sort by chromosome and position
    sort -k1,1 -k2,2n ${genotype}_all_peaks.bed > ${genotype}_all_peaks.sorted.bed
    
    # Merge peaks within ${params.merge_distance}bp
    bedtools merge -i ${genotype}_all_peaks.sorted.bed -d ${params.merge_distance} > ${genotype}_unified_regions.bed
    
    REGION_COUNT=\$(wc -l < ${genotype}_unified_regions.bed)
    echo "✅ ${genotype}: \$REGION_COUNT unified regions"
    """
}

process createConsensus {
    publishDir "${params.output_dir}", mode: 'copy'
    
    input:
    path(all_unified_beds)
    
    output:
    path("consensus_binding_regions.bed")
    path("consensus_stats.txt")
    
    script:
    """
    # Find regions bound in ALL genotypes
    # Start with first genotype
    cp ${all_unified_beds[0]} temp.bed
    
    # Intersect with each additional genotype
    for bed in ${all_unified_beds.drop(1).join(' ')}; do
        bedtools intersect -a temp.bed -b \$bed -u > temp2.bed
        mv temp2.bed temp.bed
    done
    
    # Final sort and merge
    sort -k1,1 -k2,2n temp.bed | bedtools merge > consensus_binding_regions.bed
    
    # Generate statistics
    echo "Consensus Binding Region Analysis" > consensus_stats.txt
    echo "=================================" >> consensus_stats.txt
    echo "" >> consensus_stats.txt
    echo "Parameters:" >> consensus_stats.txt
    echo "  Genotypes: ${genotypes.join(', ')}" >> consensus_stats.txt
    echo "  Targets merged: ${targets.join(', ')}" >> consensus_stats.txt
    echo "  Merge distance: ${params.merge_distance}bp" >> consensus_stats.txt
    echo "" >> consensus_stats.txt
    echo "Results:" >> consensus_stats.txt
    
    for bed in ${all_unified_beds}; do
        genotype=\$(basename \$bed _unified_regions.bed)
        count=\$(wc -l < \$bed)
        echo "  \$genotype unified regions: \$count" >> consensus_stats.txt
    done
    
    CONSENSUS_COUNT=\$(wc -l < consensus_binding_regions.bed)
    echo "  Consensus regions (in all genotypes): \$CONSENSUS_COUNT" >> consensus_stats.txt
    
    cat consensus_stats.txt
    echo ""
    echo "✅ Created consensus binding regions: \$CONSENSUS_COUNT regions"
    """
}

workflow {
    // Load overlap BED files for specified genotypes and targets
    Channel
        .fromPath("${params.overlap_dir}/*_R1_vs_R2_overlap.bed")
        .map { file ->
            def filename = file.name
            def parts = filename.replace('_R1_vs_R2_overlap.bed', '').tokenize('_')
            def genotype_val = parts[0]
            def target_val = parts[1]
            [genotype_val, target_val, file]
        }
        .filter { genotype_val, target_val, file -> 
            genotype_val in genotypes && target_val in targets 
        }
        .set { overlap_files }
    
    // Create reproducible peaks from overlaps
    reproducible = createReproduciblePeaks(overlap_files)
    
    // Group by genotype to merge targets
    reproducible
        .map { genotype_val, target_val, bed -> [genotype_val, bed] }
        .groupTuple()
        .set { grouped_by_genotype }
    
    // Merge targets per genotype
    unified = mergeTargets(grouped_by_genotype)
    
    // Collect all unified beds for consensus
    all_unified = unified.map { it[1] }.collect()
    
    // Create consensus across all genotypes
    createConsensus(all_unified)
}
