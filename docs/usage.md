# Usage Guide

## Basic Usage
```bash
nextflow run main.nf --sample_list samples.txt
```

## Advanced Options

- `--bam_dir`: Directory containing BAM files (default: /home/ec2-user/cutnrun/full_run/bams)
- `--outdir`: Output directory (default: ./results)
- `--sample_list`: File with sample names (default: samples.txt)

## Running Individual Modules

### BAM Processing Only
```bash
nextflow run modules/bam_processing.nf --sample_id SAMPLE --bam_file path/to/file.bam
```
