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
- Required for replicate averaging

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
└── averaged_bedgraphs/
    └── {genotype}_{antibody}_averaged.bedgraph
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

## Troubleshooting

### Files Going to `null/` Directory
- **Cause:** Missing `params.outdir` in module
- **Fix:** Ensure each module has `params.outdir = 'results'` at top

### Pipeline Not Resuming
- **Cause:** Changed process names
- **Fix:** Process name changes break cache; re-run without `-resume`

### Normalization Not Running
- **Cause:** Process outputs not emitted as tuples
- **Fix:** Ensure `emit: bedgraph = ...` in BAM_PROCESSING workflow

## Performance Tips

1. **Use `-resume`** when iterating on modules
2. **Don't delete `work/`** unless necessary
3. **Adjust `maxForks`** in `nextflow.config` based on your CPU count
4. **Use `stageInMode = 'symlink'`** for large BAM files (already configured)

## Running Individual Modules

### BAM Processing Only
Not currently supported - use full pipeline with `-resume` instead.
