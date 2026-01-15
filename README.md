# CUT&RUN Analysis Pipeline

Modular Nextflow pipeline for CUT&RUN data analysis from BAM files to peaks.

## Pipeline Overview
```
BAMs → [1] Process → [2] QC Stats → [3] Normalize → [4] Call Peaks
       ✓ Automated   🚧 Coming soon   🚧 Coming soon   🚧 Coming soon
```

Currently implemented: **BAM Processing Module**

## Quick Start
```bash
# Run BAM processing
nextflow run main.nf --sample_list samples.txt

# Run in background
nextflow run main.nf -bg > pipeline.log 2>&1

# Monitor progress
tail -f pipeline.log
```

## Modules

### 1. BAM Processing (✓ Complete)
- Quality filtering (MAPQ ≥ 30)
- Blacklist removal (ENCODE)
- Fragment extraction (< 1000 bp)
- Coverage track generation

**Coming soon:**
- Spike-in QC
- Peak calling (SEACR/MACS2)
- Differential binding analysis

## Input

Text file with one sample name per line (without .bam extension):
```
sample1
sample2
sample3
```

BAM files should be in the directory specified by `--bam_dir` (default: `/home/ec2-user/cutnrun/full_run/bams`)

## Output

Results in `results/bam_processing/`:
- `{sample}_final.bam` - Processed BAM file
- `{sample}_final.bam.bai` - BAM index
- `{sample}_final.clean.bedpe` - Filtered paired-end fragments
- `{sample}_final.fragments.bed` - Fragment coordinates
- `{sample}_final.fragments.sorted.bedgraph` - Coverage track
- `{sample}_run.log` - Processing log

## Performance

### Benchmarking

Tested on AWS EC2 m7i.16xlarge (64 CPUs, 247 GB RAM):

| Samples | Parallel Jobs | Wall Time | Time per Sample |
|---------|---------------|-----------|-----------------|
| 2       | 2             | 13m 25s   | ~6.5 min        |
| 72      | 8             | ~60 min*  | ~6.5 min        |

*Estimated based on test run

### Resource Usage
- CPU: 8 cores per sample
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
- `--outdir`: Output directory (default: `./results`)

## Documentation

- [Usage Guide](docs/usage.md)
- [Spike-in Normalization](docs/spike_in_normalization.md) (coming soon)

## Author

Garrett Cooper  
Emory University

## Citation

If you use this pipeline, please cite:  
Cooper et al. (2026) *Nature Communications* (in review)
