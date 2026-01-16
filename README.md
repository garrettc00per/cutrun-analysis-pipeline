# CUT&RUN Analysis Pipeline

Modular Nextflow pipeline for CUT&RUN data analysis from BAM files through peak calling and consensus analysis.

## Pipeline Overview
```
BAMs → [1] Process → [2] Normalize → [3] BigWigs → [4] Call Peaks → [5] Consensus
       ✓ Complete    ✓ Complete      ✓ Complete   ✓ Complete      ✓ Complete
```

**Implemented:**
- ✅ BAM Processing Module
- ✅ Bedgraph Normalization Module
- ✅ BigWig Generation Module
- ✅ Peak Calling Module (SEACR with adaptive filtering)
- ✅ Consensus Binding Analysis (optional workflow)

## Quick Start
```bash
# Run full pipeline
nextflow run main.nf --sample_list samples.txt

# Run in background
nextflow run main.nf --sample_list samples.txt -bg > pipeline.log 2>&1

# Monitor progress
tail -f pipeline.log

# Optional: Consensus analysis
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
- Adaptive filtering (hybrid SNR + fixed threshold)
- Replicate overlap analysis
- Comprehensive QC metrics

### 5. Consensus Analysis (✅ Complete)
- Reproducible peaks from replicate overlaps
- Merge multiple targets (e.g., protein complex subunits)
- Generate consensus regions across genotypes
- Fully customizable via command line

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
Pre-processed alignment files:
- Aligned to reference genome (e.g., GRCh38)
- Adapter-trimmed and PCR duplicates removed
- Indexed (.bam.bai files present)
- Located in directory specified by `--bam_dir`

**Naming convention:** `{genotype}_{target}_{replicate}`
- **Critical:** IgG controls required for each genotype/replicate

### normalization_factors.tsv (optional)
Tab-separated file with spike-in normalization factors:
```
Prefix	Factor
WT_IgG_R1	0.520284897
WT_H3K4me3_R1	2.429385008
...
```

## Output

### results/bam_processing/
- `{sample}_final.bam` - Processed BAM file
- `{sample}_final.fragments.bed` - Fragment coordinates
- `{sample}_final.fragments.sorted.bedgraph` - Coverage track

### results/normalized_bedgraphs/
- `{sample}_normalized.bedgraph` - Spike-in normalized coverage

### results/bigwigs_scaled/
- `{sample}_scaled.bw` - Normalized bigWig files

### results/bigwigs_igg_subtracted/
- `{sample}_IgGsubtracted.bw` - Background-subtracted bigWigs

### results/seacr_peaks/
- `{sample}_SEACR_peaks.stringent.bed` - Raw peak calls

### results/filtered_peaks/
- `{sample}_filtered.bed` - Filtered peaks (adaptive threshold)
- `{sample}_filter_summary.txt` - Filtering statistics

### results/replicate_overlaps/
- `{genotype}_{target}_R1_vs_R2_overlap.bed` - Overlapping peaks
- `replicate_overlap_summary.txt` - Summary table

### results/complex_binding/ (from consensus workflow)
- `consensus_binding_regions.bed` - Final consensus peaks
- `consensus_stats.txt` - Analysis summary

## Performance

### Benchmarking

Tested on AWS EC2 m7i.16xlarge (64 CPUs, 247 GB RAM):

| Samples | Module           | Wall Time   | Notes                    |
|---------|------------------|-------------|--------------------------|
| 2       | BAM Processing   | 13m 25s     | ~6.5 min per sample      |
| 2       | Normalization    | <1 min      | Very fast                |
| 2       | BigWig Gen       | ~2.5 min    | Fast                     |
| 60      | Peak Calling     | ~45 min     | SEACR + filtering        |
| 72      | **Full Pipeline**| **~5 hours**| End-to-end               |

### Resource Usage
- CPU: 8 cores (BAM), 4 cores (BigWig), 2 cores (peaks), 1 core (norm)
- Memory: 4-8 GB per sample
- Disk: ~2× input BAM size for temporary work files
- Parallelization: Up to 18 samples simultaneously (configurable)

## Requirements

- Nextflow (≥ 21.10.0)
- samtools (≥ 1.15)
- bedtools (≥ 2.30)
- deepTools (≥ 3.5.0)
- SEACR (v1.3) - included
- ENCODE blacklist file - included
- Genome chromosome sizes file - included

## Parameters

### Main Pipeline
- `--sample_list`: Sample names file (default: `samples.txt`)
- `--bam_dir`: BAM directory (default: `/home/ec2-user/cutnrun/full_run/bams`)
- `--norm_factors`: Normalization factors (default: `normalization_factors.tsv`)
- `--outdir`: Output directory (default: `results`)
- `--filter_method`: Peak filtering strategy (default: `hybrid`)
- `--fixed_threshold`: Minimum peak intensity (default: `5.0`)
- `--test_mode`: Test with WT only (default: `false`)

### Consensus Analysis
- `--genotypes`: Comma-separated genotypes (default: `WT,I315I`)
- `--targets_to_merge`: Comma-separated targets (default: `SMARCA4,SMARCB1,SMARCE1`)
- `--merge_distance`: Peak merge distance in bp (default: `100`)

## Key Features

**Adaptive Peak Filtering**
- Calculates IgG background statistics per sample
- Hybrid threshold: max(SNR, fixed threshold)
- Handles biological negatives (validated on nonsense mutants)

**Modular Design**
- Independent, reusable modules
- Optional analysis workflows
- Easy to customize and extend

**Production Ready**
- Comprehensive error handling
- Detailed QC reports
- Resume capability for failed runs
- Resource-aware parallelization

## Iterative Development

The pipeline supports Nextflow's `-resume` flag:
```bash
# First run
nextflow run main.nf --sample_list samples.txt

# Modify a module
nano modules/peak_calling.nf

# Re-run - only changed parts execute
nextflow run main.nf --sample_list samples.txt -resume
```

## Documentation

- [Usage Guide](docs/usage.md) - Detailed examples and troubleshooting
- [Development Notes](DEVELOPMENT.md) - Technical details and design decisions

## Author

Garrett Cooper  
Emory University  
PhD, Genetics and Molecular Biology

## Citation

If you use this pipeline, please cite:  
Cooper et al. (2026) *Nature Communications* (in review)
