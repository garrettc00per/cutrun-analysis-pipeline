# Usage Guide

## Basic Usage

### Full Pipeline (BAMs → Peaks)
```bash
nextflow run main.nf
```

### Test Mode (WT samples only)
```bash
nextflow run main.nf --test_mode true
```

### Background Execution
```bash
nextflow run main.nf -bg > pipeline.log 2>&1
tail -f pipeline.log  # Monitor progress
```

### Consensus Binding Analysis (Optional)
```bash
# After main pipeline completes
nextflow run workflows/complex_binding.nf \
  --genotypes WT,I315I \
  --targets_to_merge SMARCA4,SMARCB1,SMARCE1
```

## Input Requirements

### Main Pipeline

Normalized bedgraph files in `results/normalized_bedgraphs/` with naming pattern:
```
{genotype}_{target}_{replicate}_normalized.bedgraph
```

Examples:
```
WT_SMARCB1_R1_normalized.bedgraph
WT_SMARCB1_R2_normalized.bedgraph
WT_IgG_R1_normalized.bedgraph
WT_IgG_R2_normalized.bedgraph
I315I_H3K27me3_R1_normalized.bedgraph
I315I_H3K27me3_R2_normalized.bedgraph
```

**Naming convention:** `{genotype}_{target}_{replicate}`
- **Genotype**: WT, I315I, W281P, etc.
- **Target**: SMARCB1, H3K27me3, IgG, etc.
- **Replicate**: R1, R2

**Important:** IgG controls must be present for each genotype/replicate combination.

### Consensus Analysis

Replicate overlap BED files from main pipeline in `results/replicate_overlaps/`

## Parameters

### Main Pipeline Parameters

| Parameter           | Default                                          | Description                           |
|---------------------|--------------------------------------------------|---------------------------------------|
| `--bedgraph_dir`    | `${projectDir}/results/normalized_bedgraphs`     | Input bedgraph directory              |
| `--outdir`          | `results`                                        | Output directory                      |
| `--seacr_script`    | `${projectDir}/SEACR_1.3.sh`                     | Path to SEACR script                  |
| `--test_mode`       | `false`                                          | Run with WT samples only              |
| `--filter_method`   | `hybrid`                                         | Peak filtering strategy               |
| `--fixed_threshold` | `5.0`                                            | Minimum peak intensity                |
| `--snr_multiplier`  | `2.0`                                            | SNR threshold multiplier              |

### Consensus Analysis Parameters

| Parameter            | Default                           | Description                          |
|----------------------|-----------------------------------|--------------------------------------|
| `--overlap_dir`      | `results/replicate_overlaps`      | Input overlap BED directory          |
| `--output_dir`       | `results/complex_binding`         | Output directory                     |
| `--genotypes`        | `WT,I315I`                        | Comma-separated genotypes            |
| `--targets_to_merge` | `SMARCA4,SMARCB1,SMARCE1`         | Comma-separated targets              |
| `--merge_distance`   | `100`                             | Distance for merging peaks (bp)      |

## Peak Filtering Strategies

The pipeline offers three filtering methods:

### 1. Hybrid (Recommended - Default)
```bash
nextflow run main.nf --filter_method hybrid
```
Uses `max(SNR threshold, fixed threshold)`:
- Adapts to IgG quality (uses higher SNR if IgG is noisy)
- Protects against biological negatives (enforces minimum threshold)
- Best for diverse datasets with varying IgG quality

### 2. Signal-to-Noise Ratio (SNR) Only
```bash
nextflow run main.nf --filter_method snr --snr_multiplier 3.0
```
Uses only `mean(IgG) + N×SD(IgG)`:
- Purely statistical approach
- Adapts to each sample's IgG background
- May be too lenient for clean IgG controls

### 3. Fixed Threshold Only
```bash
nextflow run main.nf --filter_method fixed --fixed_threshold 10.0
```
Uses only a fixed intensity cutoff:
- Simple, consistent threshold
- Doesn't account for IgG quality
- Good for well-controlled experiments

## Output Structure
```
results/
├── seacr_peaks/
│   └── {sample}_SEACR_peaks.stringent.bed
├── filtered_peaks/
│   ├── {sample}_filtered.bed
│   └── {sample}_filter_summary.txt
├── peak_summary/
│   └── peak_calling_summary.txt
├── replicate_overlaps/
│   ├── {genotype}_{target}_R1_vs_R2_overlap.bed
│   ├── {genotype}_{target}_overlap_stats.txt
│   └── replicate_overlap_summary.txt
└── complex_binding/                          # From consensus analysis
    ├── reproducible_peaks/
    │   └── {genotype}_{target}_reproducible_peaks.bed
    ├── unified_regions/
    │   └── {genotype}_unified_regions.bed
    ├── consensus_binding_regions.bed
    └── consensus_stats.txt
```

## Module Outputs

### Peak Calling Module

**SEACR Peaks** (`seacr_peaks/`)
- Raw peak calls using stringent mode with IgG controls
- Format: BED file with peak coordinates and max intensity

**Filtered Peaks** (`filtered_peaks/`)
- Peaks passing hybrid/SNR/fixed threshold
- Per-sample filtering summaries showing original vs filtered counts

**Peak Summary** (`peak_summary/`)
- Complete table of all samples with:
  - Original peak count
  - Filtered peak count
  - Percent retained
  - SNR threshold
  - Fixed threshold
  - Applied threshold
  - Filter method used

**Replicate Overlaps** (`replicate_overlaps/`)
- Peaks found in both R1 and R2
- Overlap statistics (counts and percentages)
- Summary table for all conditions

### Consensus Binding Analysis

**Reproducible Peaks** (`reproducible_peaks/`)
- Union of overlapping R1 and R2 peaks per condition
- High-confidence peaks for downstream analysis

**Unified Regions** (`unified_regions/`)
- Merged peaks across specified targets (e.g., SWI/SNF subunits)
- Represents regions bound by any of the merged targets

**Consensus Regions** (main output)
- Final consensus peaks found across all specified genotypes
- Represents high-confidence binding sites common to all conditions

## Usage Examples

### Example 1: Standard Run (All Samples)
```bash
# Run peak calling on all samples
nextflow run main.nf

# Check results
cat results/peak_summary/peak_calling_summary.txt
cat results/replicate_overlaps/replicate_overlap_summary.txt
```

### Example 2: Test with WT Only
```bash
# Quick test with just WT samples
nextflow run main.nf --test_mode true

# Verify outputs
ls results/filtered_peaks/WT_*
```

### Example 3: Custom Filtering
```bash
# Use more stringent SNR filtering
nextflow run main.nf \
  --filter_method snr \
  --snr_multiplier 3.0

# Or use higher fixed threshold
nextflow run main.nf \
  --filter_method fixed \
  --fixed_threshold 10.0
```

### Example 4: SWI/SNF Complex Analysis
```bash
# First run main pipeline
nextflow run main.nf

# Then analyze SWI/SNF complex binding
nextflow run workflows/complex_binding.nf \
  --genotypes WT,I315I \
  --targets_to_merge SMARCA4,SMARCB1,SMARCE1

# Check output
wc -l results/complex_binding/consensus_binding_regions.bed
cat results/complex_binding/consensus_stats.txt
```

### Example 5: Compare Mutants
```bash
# Analyze SMARCB1 binding in different mutants
nextflow run workflows/complex_binding.nf \
  --genotypes W281P,W281X,I315X \
  --targets_to_merge SMARCB1 \
  --output_dir results/smarcb1_mutant_comparison
```

### Example 6: Histone Mark Co-localization
```bash
# Find regions with both H3K27me3 and H3K4me3
nextflow run workflows/complex_binding.nf \
  --genotypes WT,I315I \
  --targets_to_merge H3K27me3,H3K4me3 \
  --merge_distance 500 \
  --output_dir results/bivalent_chromatin
```

### Example 7: Single Subunit Analysis
```bash
# Just analyze SMARCA4 across two genotypes
nextflow run workflows/complex_binding.nf \
  --genotypes WT,I315R \
  --targets_to_merge SMARCA4 \
  --output_dir results/smarca4_wt_vs_i315r
```

## Advanced Usage

### Resume After Failure
```bash
nextflow run main.nf -resume
```
Only re-runs failed or incomplete processes.

### Clean Up Between Runs
```bash
# Remove work files and start fresh
rm -rf work/ .nextflow*

# Remove specific output directories
rm -rf results/seacr_peaks results/filtered_peaks
```

### Check Pipeline Progress
```bash
# View real-time logs
tail -f .nextflow.log

# Monitor resource usage
htop  # or top

# Check completed processes
nextflow log
```

### Generate Execution Reports
```bash
nextflow run main.nf \
  -with-report report.html \
  -with-timeline timeline.html \
  -with-dag flowchart.png
```

## Quality Control

### Inspect IgG Statistics
```bash
# Check IgG background levels
grep "IgG stats" .nextflow.log

# View individual sample thresholds
cat results/filtered_peaks/*_filter_summary.txt
```

### Assess Replicate Concordance
```bash
# View overlap summary
column -t results/replicate_overlaps/replicate_overlap_summary.txt

# Flag poor replicates (< 50% overlap)
awk -F'\t' 'NR>1 && ($6<50 || $7<50)' \
  results/replicate_overlaps/replicate_overlap_summary.txt
```

### Check Peak Counts
```bash
# Count peaks per sample
wc -l results/filtered_peaks/*.bed

# Identify outliers (very few or very many peaks)
for f in results/filtered_peaks/*_filtered.bed; do
  echo "$(wc -l < "$f") $f"
done | sort -n
```

## Troubleshooting

### No Peaks Called
**Symptoms:** Empty or very few peaks in output files

**Possible causes:**
1. IgG controls missing or incorrectly named
2. Threshold too stringent
3. Poor data quality

**Solutions:**
```bash
# Check IgG files exist
ls results/normalized_bedgraphs/*IgG*

# Try more lenient filtering
nextflow run main.nf --filter_method fixed --fixed_threshold 2.0

# Check individual sample quality
cat results/filtered_peaks/{sample}_filter_summary.txt
```

### Poor Replicate Overlap
**Symptoms:** Low overlap percentages in `replicate_overlap_summary.txt`

**Possible causes:**
1. Technical variation between replicates
2. Low sequencing depth
3. Threshold too stringent

**Solutions:**
- Inspect raw SEACR calls: `wc -l results/seacr_peaks/*.bed`
- Lower filtering threshold
- Check if issue is specific to certain samples/antibodies

### Consensus Analysis Finds No Peaks
**Symptoms:** Empty or very small `consensus_binding_regions.bed`

**Possible causes:**
1. Genotypes have very different binding patterns
2. Wrong targets specified
3. Merge distance too small

**Solutions:**
```bash
# Check individual genotype unified regions
wc -l results/complex_binding/unified_regions/*.bed

# Try larger merge distance
nextflow run workflows/complex_binding.nf \
  --genotypes WT,I315I \
  --targets_to_merge SMARCA4,SMARCB1,SMARCE1 \
  --merge_distance 500
```

### Pipeline Crashes or Stalls
**Possible causes:**
1. Insufficient memory
2. Missing dependencies
3. Corrupt input files

**Solutions:**
```bash
# Check system resources
free -h
df -h

# Verify SEACR is accessible
bash SEACR_1.3.sh --help

# Test with small subset
nextflow run main.nf --test_mode true
```

## Performance Tips

1. **Use `-resume`** when iterating - saves hours on large datasets
2. **Adjust `maxForks`** in config based on available CPUs
3. **Clean `work/` periodically** to save disk space
4. **Monitor memory usage** during peak calling
5. **Run test mode first** to validate settings before full run

## Downstream Analysis

### Load Peaks in IGV
```bash
# Filtered peaks are ready for IGV
# Load: results/filtered_peaks/{sample}_filtered.bed

# For cleaner visualization, use reproducible peaks
# Load: results/complex_binding/reproducible_peaks/*.bed
```

### Differential Binding Analysis
```bash
# Coming soon: DiffBind module
# For now, export filtered peaks for external tools
```

### Motif Analysis
```bash
# Extract FASTA sequences from consensus regions
bedtools getfasta -fi genome.fa \
  -bed results/complex_binding/consensus_binding_regions.bed \
  -fo consensus_sequences.fa

# Run MEME, HOMER, or other motif tools
```

### Annotation
```bash
# Annotate peaks with nearest genes
# Use ChIPseeker, HOMER, or bedtools closest
```
