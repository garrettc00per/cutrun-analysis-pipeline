# CUT&RUN Analysis Pipeline

Modular Nextflow pipeline for CUT&RUN data analysis from BAM files to peaks.

## Pipeline Overview
```
BAMs → [1] Process → [2] Normalize → [3] Call Peaks → [4] Diff Binding
       ✓ Complete    ✓ Complete      🚧 Coming soon   🚧 Coming soon
```

**Currently implemented:**
- ✅ BAM Processing Module
- ✅ Bedgraph Normalization Module

## Quick Start
```bash
# Prepare input files
# 1. samples.txt - one sample name per line
# 2. normalization_factors.tsv - spike-in normalization factors

# Run full pipeline
nextflow run main.nf --sample_list samples.txt

# Run in background
nextflow run main.nf -bg > pipeline.log 2>&1

# Monitor progress
tail -f pipeline.log
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

**Coming soon:**
- Automated spike-in alignment (BWA)
- Peak calling (SEACR)
- Differential binding analysis (DiffBind)

## Input

### samples.txt
Text file with one sample name per line (without .bam extension):
```
WT_IgG_R1
WT_H3K4me3_R1
WT_H3K27me3_R1
```

### normalization_factors.tsv
Tab-separated file with spike-in normalization factors:
```
Prefix	Factor
WT_IgG_R1	0.520284897
WT_H3K4me3_R1	2.429385008
WT_H3K27me3_R1	17.86575146
```

BAM files should be in the directory specified by `--bam_dir` (default: `/home/ec2-user/cutnrun/full_run/bams`)

## Output

### results/bam_processing/
- `{sample}_final.bam` - Processed BAM file
- `{sample}_final.bam.bai` - BAM index
- `{sample}_final.clean.bedpe` - Filtered paired-end fragments
- `{sample}_final.fragments.bed` - Fragment coordinates
- `{sample}_final.fragments.sorted.bedgraph` - Coverage track
- `{sample}_run.log` - Processing log

### results/normalized_bedgraphs/
- `{sample}_normalized.bedgraph` - Spike-in normalized coverage

### results/baseline_normalized/
- `{sample}_baseline.bedgraph` - Baseline-corrected for visualization

### results/averaged_bedgraphs/
- `{genotype}_{antibody}_averaged.bedgraph` - Averaged replicates

## Performance

### Benchmarking

Tested on AWS EC2 m7i.16xlarge (64 CPUs, 247 GB RAM):

| Samples | Module           | Wall Time   | Notes                    |
|---------|------------------|-------------|--------------------------|
| 2       | BAM Processing   | 13m 25s     | ~6.5 min per sample      |
| 2       | Normalization    | <1 min      | Very fast                |
| 2       | **Full Pipeline**| **16m 19s** | End-to-end               |
| 72      | Full Pipeline    | ~70 min*    | Estimated                |

*Estimated based on test runs

### Resource Usage
- CPU: 8 cores per sample (BAM processing), 1 core (normalization)
- Memory: Varies by BAM size
- Disk: Temporary work files require ~2x input BAM size
- Parallelization: 8 samples simultaneously (configurable via `maxForks` in config)

## Requirements

- Nextflow (≥ 21.10.0)
- samtools (≥ 1.15)
- bedtools (≥ 2.30)
- ENCODE blacklist file (included)
- Genome chromosome sizes file (included)

## Parameters

- `--sample_list`: Text file with sample names (default: `samples.txt`)
- `--bam_dir`: Directory containing input BAM files (default: `/home/ec2-user/cutnrun/full_run/bams`)
- `--norm_factors`: Normalization factors file (default: `normalization_factors.tsv`)
- `--outdir`: Output directory (default: `./results`)

## Iterative Development

The pipeline supports Nextflow's `-resume` flag for efficient iteration:
```bash
# First run
nextflow run main.nf --sample_list samples.txt

# Modify a module
nano modules/bedgraph_normalization.nf

# Re-run - only changed parts execute
nextflow run main.nf --sample_list samples.txt -resume
```

## Documentation

- [Usage Guide](docs/usage.md)
- [Development Notes](DEVELOPMENT.md)

## Author

Garrett Cooper  
Emory University

## Citation

If you use this pipeline, please cite:  
Cooper et al. (2026) *Nature Communications* (in review)
