# Development Notes

## Pipeline Architecture

### Current Structure
- **main.nf**: Orchestrator that calls modules
- **modules/bam_processing.nf**: BAM QC and fragment extraction
- Outputs go to `results/bam_processing/`

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

**Performance:**
- ~6.5 minutes per sample on AWS EC2 (64 CPUs)
- 2 samples in parallel: 13m 25s
- Full 72 samples: ~60 minutes estimated

### Known Issues / Gotchas

**bedtools warnings are normal:**
- "Query marked as paired but mate does not occur next to it" 
- This happens because filtering removes some read pairs
- Safe to ignore - pipeline works correctly

**File staging:**
- Originally was copying BAMs to work dirs (very slow)
- Fixed with `stageInMode = 'symlink'` in config
- Without this, pipeline is 5-10x slower

### Future Development Plans

**Phase 1: Spike-in QC Module** (Next)
- Collect alignment statistics (human vs E.coli)
- Output: CSV with percentages per sample
- User manually calculates normalization factors

**Phase 2: Peak Calling Module**
- SEACR for CUT&RUN data
- Optional: Apply spike-in normalization factors
- Generate peak BED files

**Phase 3: Differential Binding** (Optional)
- DiffBind or csaw
- Compare conditions

### Testing

**Quick test (2 samples):**
```bash
head -2 samples.txt > test_samples.txt
time nextflow run main.nf --sample_list test_samples.txt
```

**Full run (72 samples):**
```bash
nextflow run main.nf -bg > pipeline.log 2>&1
tail -f pipeline.log
```

**Clean up between runs:**
```bash
rm -rf work/ .nextflow* results/
```

### Dependencies

- Nextflow 25.10.2
- Java 17 (installed via conda)
- samtools, bedtools (available in PATH)
- ENCODE blacklist: ENCFF356LFX.bed.gz (8KB)
- Genome sizes: GRCh38.p13.chrom.sizes (11KB)

### Git/GitHub

**Repo:** https://github.com/garrettc00per/cutrun-analysis-pipeline  
(or https://github.com/garrett_c00per/cutrun-analysis-pipeline)

**Key files NOT committed (in .gitignore):**
- work/
- .nextflow*
- results/
- *.bam files
- test_samples.txt

### Pipeline Context

This pipeline was developed for SMARCB1 variant CUT&RUN data analysis as part of PhD dissertation work. Focus on chromatin remodeling complex binding (SMARCB1, SMARCA4, SMARCE1) and histone marks (H3K27me3, H3K4me3).

**Data location on AWS:** `/home/ec2-user/cutnrun/full_run/bams/`  
**72 samples total:** 6 conditions × 4 targets × 2 replicates + IgG controls
