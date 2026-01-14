# CUT&RUN BAM Processing Pipeline

Nextflow pipeline for processing paired-end CUT&RUN sequencing data from raw BAM files to fragment coverage tracks.

## Overview

This pipeline processes CUT&RUN BAM files through quality filtering, blacklist removal, and fragment extraction.

**Processing steps:**
1. Sort BAM by query name
2. Filter for high-quality paired reads (MAPQ ≥ 30)
3. Remove ENCODE blacklisted regions
4. Remove non-chromosomal reads (chrM, chrUn, random, chrEBV)
5. Extract fragment coordinates (< 1000 bp)
6. Generate coverage tracks (bedGraph)

## Requirements

- Nextflow (≥ 21.10.0)
- samtools (≥ 1.15)
- bedtools (≥ 2.30)

## Quick Start
```bash
# Run pipeline
nextflow run main.nf

# Resume if interrupted
nextflow run main.nf -resume

# Run in background
nextflow run main.nf -bg > pipeline.log 2>&1
```

## Input

- BAM files in `/home/ec2-user/cutnrun/full_run/bams/`
- Sample names in `samples.txt`
- ENCODE blacklist: `ENCFF356LFX.bed.gz`
- Genome sizes: `GRCh38.p13.chrom.sizes`

## Output

Results in `processed_bams_pe/`:
- `{sample}_final.bam` - Processed BAM
- `{sample}_final.bam.bai` - BAM index
- `{sample}_final.clean.bedpe` - Filtered fragments
- `{sample}_final.fragments.bed` - Fragment coordinates  
- `{sample}_final.fragments.sorted.bedgraph` - Coverage track
- `{sample}_run.log` - Processing log

## Performance

- ~2.5 minutes per sample
- 18 samples processed in parallel
- ~3 hours for 72 samples

## Author

Garrett Cooper  
Emory University
