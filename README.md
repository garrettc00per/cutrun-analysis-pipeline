# CUT&RUN Analysis Pipeline

Modular Nextflow pipeline for CUT&RUN data analysis from BAM files through peak calling and consensus peak analysis.

## Pipeline Overview
```
BAMs → [1] Process → [2] Normalize → [3] BigWigs → [4] Call Peaks → [5] Consensus Analysis
       ✓ Complete    ✓ Complete      ✓ Complete   ✓ Complete      ✓ Complete
```

**Implemented:**
- ✅ BAM Processing Module
- ✅ Bedgraph Normalization Module
- ✅ BigWig Generation Module (with IgG subtraction and replicate averaging)
- ✅ Peak Calling Module (SEACR with adaptive filtering and replicate overlap)
- ✅ Consensus Binding Analysis Workflow (optional, project-specific)

## Quick Start

### Full Pipeline (BAMs → Peaks)
```bash
# Run complete pipeline
nextflow run main.nf --sample_list samples.txt

# Run in background
nextflow run main.nf --sample_list samples.txt -bg > pipeline.log 2>&1
tail -f pipeline.log
```

### Consensus Binding Analysis (optional)
```bash
# After main pipeline completes, create consensus peaks
nextflow run workflows/complex_binding.nf \
  --genotypes WT,I315I \
  --targets_to_merge SMARCA4,SMARCB1,SMARCE1
```

## Modules

### 1. BAM Processing (✅ Complete)
- Quality filtering (MAPQ ≥ 30)
- Blacklist removal (ENCODE)
- Fragment extraction (< 1000 bp)
- Coverage track generation (bedGraph)

### 2. Bedgraph Normalization (✅ Complete)
- Spike-in normalization
- Baseline correction for visualization
- Replicate averaging (R1 + R2)

### 3. BigWig Generation (✅ Complete)
- Scaled bigWig creation with spike-in factors
- IgG background subtraction
- Replicate averaging for visualization

### 4. Peak Calling (✅ Complete)
- SEACR peak calling with IgG controls
- Adaptive filtering:
  - Signal-to-noise ratio (SNR) calculation from IgG
  - Hybrid threshold (max of SNR or fixed threshold)
  - Configurable filtering strategy
- Replicate overlap analysis
- Quality metrics and summary reports

### 5. Consensus Binding Analysis (✅ Optional Workflow)
- Create reproducible peaks from replicate overlaps
- Merge multiple targets (e.g., protein complex subunits)
- Generate consensus regions across genotypes
- Fully customizable for different biological questions

## Input

### samples.txt
Text file with one sample name per line (without .bam extension):
```
I315I_IgG_R1
I315I_IgG_R2
I315I_SMARCB1_R1
I315I_SMARCB1_R2
WT_IgG_R1
WT_IgG_R2
WT_SMARCB1_R1
WT_SMARCB1_R2
...
```

### BAM Files
Pre-processed alignment files required for each sample:
- Aligned to reference genome (e.g., GRCh38.p13)
- Adapter-trimmed
- PCR duplicates removed
- Indexed (.bam.bai files present)
- Located in directory specified by `--bam_dir`

**Sample naming convention:** `{genotype}_{target}_{replicate}`
- **Genotype**: WT, I315I, I315R, I315X, W281P, W281X
- **Target**: SMARCB1, SMARCA4, SMARCE1, H3K27me3, H3K4me3, IgG
- **Replicate**: R1, R2

**Critical:** IgG controls must be present for each genotype/replicate combination for peak calling to work.

### normalization_factors.tsv (optional)
Tab-separated file with spike-in normalization factors:
```
Prefix	Factor
WT_IgG_R1	0.520284897
WT_H3K4me3_R1	2.429385008
WT_H3K27me3_R1	17.86575146
```

If not provided, pipeline will skip spike-in normalization step.

## Output

### Main Pipeline Outputs

#### results/bam_processing/
- `{sample}_final.bam` - Processed BAM file
- `{sample}_final.bam.bai` - BAM index
- `{sample}_final.fragments.bed` - Fragment coordinates
- `{sample}_final.fragments.sorted.bedgraph` - Coverage track

#### results/normalized_bedgraphs/
- `{sample}_normalized.bedgraph` - Spike-in normalized coverage

#### results/seacr_peaks/
- `{sample}_SEACR_peaks.stringent.bed` - Raw SEACR peak calls

#### results/filtered_peaks/
- `{sample}_filtered.bed` - Filtered peaks (hybrid threshold applied)
- `{sample}_filter_summary.txt` - Filtering statistics per sample

#### results/peak_summary/
- `peak_calling_summary.txt` - Complete filtering summary for all samples

#### results/replicate_overlaps/
- `{genotype}_{target}_R1_vs_R2_overlap.bed` - Overlapping peaks between replicates
- `{genotype}_{target}_overlap_stats.txt` - Overlap statistics per condition
- `replicate_overlap_summary.txt` - Summary table of all overlaps

### Consensus Analysis Outputs

#### results/complex_binding/reproducible_peaks/
- `{genotype}_{target}_reproducible_peaks.bed` - Peaks found in both replicates

#### results/complex_binding/unified_regions/
- `{genotype}_unified_regions.bed` - Merged peaks across specified targets

#### results/complex_binding/
- `consensus_binding_regions.bed` - Final consensus peaks across genotypes
- `consensus_stats.txt` - Analysis summary with statistics

## Peak Calling Configuration

### Filtering Strategies

The pipeline uses a **hybrid filtering approach** by default:
```groovy
params {
    filter_method = 'hybrid'        // Options: 'snr', 'fixed', 'hybrid'
    fixed_threshold = 5.0           // Minimum peak intensity
    snr_multiplier = 2.0            // SD multiplier for SNR method
}
```

**Filter Methods:**
- `snr`: Use only signal-to-noise (mean + 2×SD from IgG)
- `fixed`: Use only fixed threshold (5.0 by default)
- `hybrid`: Use max(SNR threshold, fixed threshold) - **recommended**

### Why Hybrid?

The hybrid approach:
- Protects against noisy IgG controls (high SNR kicks in)
- Protects against biological negatives (fixed threshold kicks in)
- Validated on nonsense mutants with no functional protein

## Consensus Analysis Examples

### SWI/SNF Complex Binding
```bash
nextflow run workflows/complex_binding.nf \
  --genotypes WT,I315I \
  --targets_to_merge SMARCA4,SMARCB1,SMARCE1
```

### Compare Mutants for Single Protein
```bash
nextflow run workflows/complex_binding.nf \
  --genotypes W281P,W281X \
  --targets_to_merge SMARCB1 \
  --output_dir results/smarcb1_mutants
```

### Histone Mark Co-occurrence
```bash
nextflow run workflows/complex_binding.nf \
  --genotypes WT,I315I \
  --targets_to_merge H3K27me3,H3K4me3 \
  --merge_distance 500 \
  --output_dir results/histone_marks
```

### Custom Analysis
```bash
nextflow run workflows/complex_binding.nf \
  --overlap_dir results/replicate_overlaps \
  --output_dir results/my_custom_analysis \
  --genotypes WT,I315R,I315X \
  --targets_to_merge SMARCA4,SMARCE1 \
  --merge_distance 150
```

## Performance

### Benchmarking

Tested on AWS EC2 m7i.16xlarge (64 CPUs, 247 GB RAM):

| Samples | Module              | Wall Time      | Notes                    |
|---------|---------------------|----------------|--------------------------|
| 2       | BAM Processing      | 13m 25s        | ~6.5 min per sample      |
| 2       | Normalization       | <1 min         | Very fast                |
| 2       | BigWig Gen          | ~2.5 min       | Fast                     |
| 60      | Peak Calling        | ~45 min        | SEACR + filtering        |
| 72      | **Full Pipeline**   | **~5 hours**   | End-to-end               |
| -       | Consensus Analysis  | <5 min         | Post-processing          |

### Resource Usage
- **BAM Processing**: 8 CPUs per sample, 4-8 GB memory
- **Peak Calling**: 2 CPUs per sample, 4 GB memory
- **Parallelization**: Up to 18 samples simultaneously (configurable)
- **Storage**: Intermediate files in `work/` directory (automatically cleaned with `-resume`)

## Requirements

- Nextflow (≥ 21.10.0)
- samtools (≥ 1.15)
- bedtools (≥ 2.30)
- deepTools (≥ 3.5.0)
- SEACR (v1.3)
- ENCODE blacklist file (included)
- Genome chromosome sizes file (included)

## Parameters

### Main Pipeline
- `--sample_list`: List of sample names (default: `samples.txt`)
- `--bam_dir`: Directory containing BAM files (default: `/home/ec2-user/cutnrun/full_run/bams`)
- `--norm_factors`: Spike-in normalization factors (default: `normalization_factors.tsv`)
- `--outdir`: Output directory (default: `results`)
- `--seacr_script`: Path to SEACR script (default: `${projectDir}/SEACR_1.3.sh`)
- `--test_mode`: Test with WT samples only (default: `false`)
- `--filter_method`: Peak filtering strategy (default: `hybrid`)
- `--fixed_threshold`: Minimum peak intensity (default: `5.0`)
- `--snr_multiplier`: SNR threshold multiplier (default: `2.0`)

### Consensus Analysis
- `--overlap_dir`: Overlap BED directory (default: `results/replicate_overlaps`)
- `--output_dir`: Output directory (default: `results/complex_binding`)
- `--genotypes`: Comma-separated genotypes (default: `WT,I315I`)
- `--targets_to_merge`: Comma-separated targets (default: `SMARCA4,SMARCB1,SMARCE1`)
- `--merge_distance`: Distance for merging peaks (default: `100`)

## Pipeline Architecture
```
cutrun-analysis-pipeline/
├── main.nf                      # Main pipeline orchestrator
├── nextflow.config              # Configuration
├── modules/
│   ├── bam_processing.nf        # BAM QC and processing
│   ├── bedgraph_normalization.nf# Spike-in normalization
│   ├── bigwig_generation.nf     # BigWig creation
│   └── peak_calling.nf          # Peak calling and filtering
├── workflows/
│   └── complex_binding.nf       # Consensus analysis (optional)
├── SEACR_1.3.sh                 # Peak caller script
├── ENCFF356LFX.bed.gz           # ENCODE blacklist
└── GRCh38.p13.chrom.sizes       # Genome sizes
```

## Iterative Development

The pipeline supports Nextflow's `-resume` flag for efficient iteration:
```bash
# First run
nextflow run main.nf --sample_list samples.txt

# Modify a module
nano modules/peak_calling.nf

# Re-run - only changed parts execute
nextflow run main.nf --sample_list samples.txt -resume
```

## Documentation

- [Usage Guide](docs/usage.md)
- [Development Notes](DEVELOPMENT.md)

## Key Features

### Adaptive Peak Filtering
- Calculates IgG background statistics per sample
- Hybrid threshold combines statistical and biological criteria
- Handles biological negatives (e.g., nonsense mutants)
- Validated on 72 samples across 6 genotypes

### Modular Design
- Independent, reusable modules
- Optional analysis workflows
- Easy to add new modules or customize existing ones

### Production Ready
- Comprehensive error handling
- Detailed logging and reporting
- Resume capability for failed runs
- Resource-aware parallelization

## Author

Garrett Cooper  
Emory University  
PhD Candidate, Genetics and Molecular Biology

## Citation

If you use this pipeline, please cite:  
Cooper et al. (2026) *Nature Communications* (in review)
