# Usage Guide

## Basic Usage

### Full Pipeline
```bash
nextflow run main.nf --sample_list samples.txt
```

### With Custom Parameters
```bash
nextflow run main.nf \
  --sample_list my_samples.txt \
  --bam_dir /path/to/bams \
  --norm_factors my_normalization_factors.tsv \
  --outdir my_results
```

### Background Execution
```bash
nextflow run main.nf -bg > pipeline.log 2>&1
tail -f pipeline.log  # Monitor progress
```

## Input File Formats

### samples.txt
One sample name per line (no .bam extension):
```
WT_IgG_R1
WT_H3K4me3_R1
WT_H3K27me3_R1
```

### normalization_factors.tsv
Tab-separated with header:
```
Prefix	Factor
WT_IgG_R1	0.520284897
WT_H3K4me3_R1	2.429385008
```

**Sample naming convention:** `{genotype}_{antibody}_{replicate}`
- Example: `WT_H3K27me3_R1`
- Required for replicate averaging and IgG matching

## Parameters

| Parameter       | Default                                     | Description                          |
|-----------------|---------------------------------------------|--------------------------------------|
| `--sample_list` | `samples.txt`                               | List of sample names                 |
| `--bam_dir`     | `/home/ec2-user/cutnrun/full_run/bams`      | Directory with BAM files             |
| `--norm_factors`| `normalization_factors.tsv`                 | Spike-in normalization factors       |
| `--outdir`      | `results`                                   | Output directory                     |

## Output Structure
```
results/
├── bam_processing/
│   ├── {sample}_final.bam
│   ├── {sample}_final.bam.bai
│   ├── {sample}_final.clean.bedpe
│   ├── {sample}_final.fragments.bed
│   ├── {sample}_final.fragments.sorted.bedgraph
│   └── {sample}_run.log
├── normalized_bedgraphs/
│   └── {sample}_normalized.bedgraph
├── baseline_normalized/
│   └── {sample}_baseline.bedgraph
├── averaged_bedgraphs/
│   └── {genotype}_{antibody}_averaged.bedgraph
├── bigwigs_scaled/
│   └── {sample}_scaled.bw
├── bigwigs_igg_subtracted/
│   └── {sample}_IgGsubtracted.bw
└── bigwigs_averaged/
    └── {genotype}_{antibody}_avg50bp.bw
```

## Advanced Usage

### Resume After Failure
```bash
nextflow run main.nf --sample_list samples.txt -resume
```

### Test with Subset
```bash
head -2 samples.txt > test_samples.txt
nextflow run main.nf --sample_list test_samples.txt
```

### Clean Up
```bash
# Remove cache and results
rm -rf work/ .nextflow* results/

# Remove only bad outputs (keep cache)
rm -rf null/
```

## Module Outputs

### BAM Processing
- Filtered, sorted BAMs ready for downstream analysis
- Fragment BED files for peak calling
- BedGraph coverage tracks

### Bedgraph Normalization
- Spike-in normalized bedGraphs for quantitative comparisons
- Baseline-corrected bedGraphs for visualization
- Averaged replicate bedGraphs (R1 + R2)

### BigWig Generation
- Scaled bigWigs for genome browser visualization
- IgG-subtracted bigWigs showing true enrichment
- Averaged bigWigs combining biological replicates

## Troubleshooting

### Files Going to `null/` Directory
- **Cause:** Missing `params.outdir` in module
- **Fix:** Ensure each module has `params.outdir = 'results'` at top

### Pipeline Not Resuming
- **Cause:** Changed process names
- **Fix:** Process name changes break cache; re-run without `-resume`

### BigWig Generation Failing
- **Cause:** deepTools not installed
- **Fix:** `conda install -c bioconda deeptools`

### IgG Subtraction Missing Samples
- **Cause:** Sample naming doesn't match expected format
- **Fix:** Ensure samples follow `{genotype}_{antibody}_{rep}` format

### Averaged BigWigs Incomplete
- **Cause:** Missing replicates or IgG subtraction incomplete
- **Fix:** Verify both R1 and R2 completed IgG subtraction successfully

## Performance Tips

1. **Use `-resume`** when iterating on modules
2. **Don't delete `work/`** unless necessary
3. **Adjust `maxForks`** in `nextflow.config` based on your CPU count
4. **Use `stageInMode = 'symlink'`** for large BAM files (already configured)
5. **BigWig generation is memory-intensive** - monitor with `top` or `htop`

## Visualization Workflows

### Load in IGV
```bash
# BigWigs can be loaded directly into IGV
# Use bigwigs_igg_subtracted/ for cleanest visualization
# Use bigwigs_averaged/ for comparing conditions
```

### Generate Plots
```bash
# Coming soon: Built-in plotting module
# For now, use deepTools or custom R/Python scripts
```
