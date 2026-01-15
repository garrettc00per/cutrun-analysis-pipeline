# Development Notes

## Pipeline Architecture

### Current Structure
- **main.nf**: Orchestrator that calls modules
- **modules/bam_processing.nf**: BAM QC and fragment extraction
- **modules/bedgraph_normalization.nf**: Spike-in normalization and replicate averaging ✅
- Outputs go to `results/bam_processing/`, `results/normalized_bedgraphs/`, `results/baseline_normalized/`, `results/averaged_bedgraphs/`

### Key Technical Decisions

**Nextflow Configuration:**
- `maxForks = 18`: Run 18 samples in parallel
- `cpus = 8`: Each sample uses 8 CPU cores
- With 64 CPUs, actually runs 8 samples at a time (8 × 8 = 64)
- `stageInMode = 'symlink'`: Don't copy BAM files, use symlinks (saves time)

**BAM Processing Steps:**
1. Sort by query name (for paired-end processing)
2. Filter: MAPQ ≥ 30, proper pairs only
3. Remove ENCODE blacklist regions (ENCFF356LFX.bed.gz)
4. Remove chrM, chrUn, random, chrEBV
5. Extract fragments < 1000 bp
6. Generate bedGraph coverage

**Normalization Steps:**
1. Apply spike-in normalization factors from `normalization_factors.tsv`
2. Baseline correction (subtract minimum signal for visualization)
3. Average biological replicates (R1 + R2)

**Performance:**
- BAM Processing: ~6.5 minutes per sample on AWS EC2 (64 CPUs)
- Normalization: <1 minute for all samples (very fast)
- 2 samples full pipeline: ~16 minutes total
- Full 72 samples: ~70 minutes estimated

### Known Issues / Gotchas

**bedtools warnings are normal:**
- "Query marked as paired but mate does not occur next to it" 
- This happens because filtering removes some read pairs
- Safe to ignore - pipeline works correctly

**File staging:**
- Originally was copying BAMs to work dirs (very slow)
- Fixed with `stageInMode = 'symlink'` in config
- Without this, pipeline is 5-10x slower

**Resume behavior:**
- Use `-resume` to skip completed steps when iterating
- Changing process names breaks cache (Nextflow sees it as new process)
- Don't delete `work/` directory unless you want to re-run everything

**publishDir requires params.outdir:**
- Each module needs `params.outdir = 'results'` at the top
- Without this, files go to `null/` directory

### Module Development Notes

**Process vs Workflow naming:**
- Process: Describes what it does (e.g., `FILTER_BAMS`)
- Workflow: Module name for import (e.g., `BAM_PROCESSING`)
- Use `as` in main.nf for clean naming: `include { RUN_BAM_PROCESSING as BAM_PROCESSING }`

**Emitting outputs:**
- Always emit as tuples: `tuple val(sample_id), path(file)`
- Normalization needs `[sample_id, file]` format to join with factors
- Use `emit:` block at end of workflow to pass data between modules

### Future Development Plans

**Phase 1: Bedgraph Normalization** ✅ COMPLETE
- Apply spike-in normalization factors
- Baseline correction for visualization
- Average biological replicates

**Phase 2: Spike-in Alignment Module** (Next)
- Align FASTQs to human + E.coli genomes with BWA
- Calculate normalization factors automatically
- Remove manual normalization_factors.tsv requirement

**Phase 3: Peak Calling Module**
- SEACR for CUT&RUN data
- Use normalized bedgraphs as input
- Generate peak BED files

**Phase 4: Differential Binding** (Optional)
- DiffBind or csaw
- Compare conditions

### Testing

**Quick test (2 samples):**
```bash
head -2 samples.txt > test_samples.txt
time nextflow run main.nf --sample_list test_samples.txt
```

**Resume after changes:**
```bash
# Modify a module
nano modules/bedgraph_normalization.nf

# Re-run - only changed parts execute
nextflow run main.nf --sample_list test_samples.txt -resume
```

**Full run (72 samples):**
```bash
nextflow run main.nf -bg > pipeline.log 2>&1
tail -f pipeline.log
```

**Clean up between runs:**
```bash
# Nuclear option - deletes everything
rm -rf work/ .nextflow* results/

# Safer - keep cache, remove bad outputs
rm -rf null/
```

### Dependencies

- Nextflow 25.10.2
- Java 17 (installed via conda)
- samtools, bedtools (available in PATH)
- ENCODE blacklist: ENCFF356LFX.bed.gz (8KB)
- Genome sizes: GRCh38.p13.chrom.sizes (11KB)
- normalization_factors.tsv (user-provided, tab-separated)

### Git/GitHub

**Repo:** https://github.com/garrettc00per/cutrun-analysis-pipeline  
(or https://github.com/garrett_c00per/cutrun-analysis-pipeline)

**Key files NOT committed (in .gitignore):**
- work/
- .nextflow*
- results/
- *.bam files
- test_samples.txt
- null/

### Pipeline Context

This pipeline was developed for SMARCB1 variant CUT&RUN data analysis as part of PhD dissertation work. Focus on chromatin remodeling complex binding (SMARCB1, SMARCA4, SMARCE1) and histone marks (H3K27me3, H3K4me3).

**Data location on AWS:** `/home/ec2-user/cutnrun/full_run/bams/`  
**72 samples total:** 6 conditions × 4 targets × 2 replicates + IgG controls
