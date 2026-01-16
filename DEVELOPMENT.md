# Development Notes

## Pipeline Architecture

### Current Structure
- **main.nf**: Main pipeline orchestrator
- **modules/bam_processing.nf**: BAM QC and fragment extraction
- **modules/bedgraph_normalization.nf**: Spike-in normalization and replicate averaging
- **modules/bigwig_generation.nf**: BigWig creation with IgG subtraction
- **modules/peak_calling.nf**: SEACR peak calling with adaptive filtering ✅
- **workflows/complex_binding.nf**: Consensus peak analysis (optional, project-specific) ✅

### Output Directories
- `results/bam_processing/` - Processed BAMs and bedGraphs
- `results/normalized_bedgraphs/` - Spike-in normalized bedGraphs
- `results/baseline_normalized/` - Visualization-ready bedGraphs
- `results/averaged_bedgraphs/` - Replicate-averaged bedGraphs
- `results/bigwigs_scaled/` - Normalized bigWigs
- `results/bigwigs_igg_subtracted/` - Background-subtracted bigWigs
- `results/bigwigs_averaged/` - Replicate-averaged bigWigs
- `results/seacr_peaks/` - Raw SEACR peak calls
- `results/filtered_peaks/` - Filtered peaks with QC stats
- `results/peak_summary/` - Peak calling summary tables
- `results/replicate_overlaps/` - R1 vs R2 overlap analysis
- `results/complex_binding/` - Consensus binding regions (from optional workflow)

## Key Technical Decisions

### Nextflow Configuration
- `maxForks = 18`: Run 18 samples in parallel
- Variable CPUs per process: 8 (BAM), 4 (BigWig), 2 (SEACR), 1 (normalization)
- With 64 CPUs, runs 8 BAM samples simultaneously (8 × 8 = 64)
- `stageInMode = 'symlink'`: Use symlinks instead of copying (critical for large BAMs)

### BAM Processing Steps
1. Sort by query name (for paired-end processing)
2. Filter: MAPQ ≥ 30, proper pairs only
3. Remove ENCODE blacklist regions (ENCFF356LFX.bed.gz)
4. Remove chrM, chrUn, random, chrEBV
5. Extract fragments < 1000 bp
6. Generate bedGraph coverage

### Peak Calling Design

**Adaptive Filtering Strategy:**
- **Hybrid threshold (default)**: `max(SNR_threshold, fixed_threshold)`
  - SNR = mean(IgG) + 2×SD(IgG)
  - Fixed = 5.0 (configurable)
- **Why hybrid?** Handles both noisy IgG controls AND biological negatives
- **Validated on nonsense mutants**: I315X/W281X SMARCB1 had no functional protein → correctly filtered

**IgG-Based Normalization:**
- Each genotype/replicate pair has matched IgG control
- IgG statistics calculated per sample
- Threshold adapts to IgG quality automatically

**Replicate Overlap:**
- Uses `bedtools intersect` to find R1 vs R2 overlaps
- Reports overlap counts and percentages
- Identifies poor concordance for QC

### Performance Benchmarks

Tested on AWS EC2 m7i.16xlarge (64 CPUs, 247 GB RAM):

| Module              | Samples | Wall Time    | Notes                          |
|---------------------|---------|--------------|--------------------------------|
| BAM Processing      | 2       | 13m 25s      | ~6.5 min per sample            |
| Normalization       | 72      | <1 min       | Very fast                      |
| BigWig Generation   | 2       | ~2.5 min     | Fast                           |
| Peak Calling        | 60      | ~45 min      | SEACR + filtering              |
| Replicate Overlap   | 30      | ~5 min       | Part of peak calling           |
| **Full Pipeline**   | **72**  | **~5 hours** | BAMs → Filtered peaks          |
| Consensus Analysis  | 6       | <5 min       | Post-processing                |

### Module Development Lessons

**Process vs Workflow Naming:**
- Process: Verb describing action (e.g., `calculateIgGStats`)
- Workflow: Module name for import (e.g., `PEAK_CALLING`)
- Workflows emit multiple outputs using `emit:` block

**Channel Operations:**
- Use `.map` to transform tuples: `[a, b, c] → [a, c]`
- Use `.filter` to subset: `{ it[0] == 'WT' }`
- Use `.join` to match on keys: `r1_peaks.join(r2_peaks, by: [0,1])`
- Use `.groupTuple` to collect by key: `[genotype, bed1], [genotype, bed2] → [genotype, [bed1, bed2]]`
- Use `.combine` for Cartesian joins with `by:` for key matching

**Parameter Handling:**
- Comma-separated CLI params need parsing: `params.genotypes.tokenize(',')`
- Check if string vs list: `params.x instanceof String ? params.x.tokenize(',') : params.x`
- Always validate required params before workflow runs

**PublishDir Best Practices:**
- Use unique parameter names to avoid conflicts (`consensus_outdir` not `outdir`)
- Always use `mode: 'copy'` for final outputs
- Tag processes for clear logging: `tag "${genotype}_${target}"`

**File Path Management:**
- `${projectDir}` for files in pipeline repo
- Absolute paths for external data
- Relative paths resolve from where `nextflow run` is executed

## Known Issues / Gotchas

### bedtools Warnings
- "Query marked as paired but mate does not occur next to it"
- Happens because filtering removes unpaired reads
- **Safe to ignore** - pipeline works correctly

### File Staging
- Originally copied BAMs to work dirs (very slow)
- Fixed with `stageInMode = 'symlink'`
- Without symlinks, pipeline is 5-10× slower

### Resume Behavior
- Use `-resume` to skip completed processes
- Changing process names breaks cache
- Don't delete `work/` unless you want to re-run everything
- Failed processes automatically retry on `-resume`

### Process Config Warnings
- "There's no process matching config selector: X"
- Occurs when running consensus workflow (doesn't have peak calling processes)
- **Safe to ignore** - just means config doesn't apply to that workflow

### Parameter Inheritance
- Workflows inherit params from main config
- Can cause conflicts if workflows use same param names
- Solution: Use unique param names (`consensus_outdir` vs `outdir`)

### Empty Consensus Results
- Happens when genotypes have very different binding patterns
- Check individual unified regions before consensus step
- Try increasing `--merge_distance` parameter

## Pipeline vs Analysis Separation

**Design principle:** Pipeline handles data processing, separate workflows handle project-specific analysis.

**Main Pipeline (general-purpose):**
- BAM processing → Normalization → BigWigs → Peak calling
- Outputs are ready for any downstream analysis
- Reusable across different CUT&RUN experiments

**Consensus Workflow (project-specific):**
- Takes peak calling outputs
- Asks specific biological questions (e.g., "where does SWI/SNF bind?")
- Fully customizable via command-line parameters
- Users can create their own analysis workflows

**Why separate?**
- Keeps pipeline clean and maintainable
- Analysis stays flexible and customizable
- Clear boundary between processing and interpretation

## Development Status

### ✅ Complete Modules
- [x] BAM Processing
- [x] Bedgraph Normalization
- [x] BigWig Generation
- [x] Peak Calling with Adaptive Filtering
- [x] Replicate Overlap Analysis
- [x] Consensus Binding Workflow

### 🚧 Future Enhancements
- [ ] Spike-in alignment module (BWA, auto-calculate normalization factors)
- [ ] Differential binding module (DiffBind or csaw)
- [ ] Peak annotation module (ChIPseeker)
- [ ] Visualization module (heatmaps, profiles)
- [ ] Quality control report (MultiQC integration)

## Testing Workflows

### Quick Test (WT samples only)
```bash
nextflow run main.nf --test_mode true
```

### Resume After Changes
```bash
# Modify a module
nano modules/peak_calling.nf

# Re-run - only changed parts execute
nextflow run main.nf -resume
```

### Full Run (72 samples)
```bash
nextflow run main.nf -bg > pipeline.log 2>&1
tail -f pipeline.log
```

### Test Consensus Workflow
```bash
# After main pipeline completes
nextflow run workflows/complex_binding.nf \
  --genotypes WT,I315I \
  --targets_to_merge SMARCA4,SMARCB1,SMARCE1
```

### Clean Up Between Runs
```bash
# Remove work files and specific outputs
rm -rf work/ .nextflow* results/seacr_peaks results/filtered_peaks

# Nuclear option - delete everything
rm -rf work/ .nextflow* results/
```

### Generate Reports
```bash
nextflow run main.nf \
  -with-report report.html \
  -with-timeline timeline.html \
  -with-dag flowchart.png
```

## Dependencies

### Required Software
- Nextflow ≥ 25.10.2
- Java 17 (for Nextflow)
- samtools ≥ 1.15
- bedtools ≥ 2.30
- deepTools ≥ 3.5.0 (for bamCoverage, bigwigCompare)
- SEACR v1.3 (included in repo)

### Reference Files (Included)
- ENCODE blacklist: `ENCFF356LFX.bed.gz` (8 KB)
- Genome sizes: `GRCh38.p13.chrom.sizes` (11 KB)

### User-Provided Files
- Normalized bedgraphs in `results/normalized_bedgraphs/`
- Naming convention: `{genotype}_{target}_{replicate}_normalized.bedgraph`
- IgG controls required for each genotype/replicate

## Git/GitHub

**Repository:** https://github.com/garrettc00per/cutrun-analysis-pipeline

### Files NOT Committed (.gitignore)
- `work/` - Nextflow work directory
- `.nextflow*` - Nextflow cache files
- `results/` - Pipeline outputs
- `*.bam`, `*.bai` - BAM files
- `*.bw` - BigWig files
- `test_samples.txt` - Test files
- `null/` - Mis-configured output directory
- `*.log` - Log files

### Commit Best Practices
- Commit after each module is working
- Tag releases: `v1.0.0`, `v1.1.0`, etc.
- Document breaking changes in commit messages
- Keep DEVELOPMENT.md and README.md in sync

## Pipeline Context

This pipeline was developed for SMARCB1 variant characterization in pediatric rhabdoid tumors as part of PhD dissertation work at Emory University.

**Research Focus:**
- Chromatin remodeling complex binding (SMARCB1, SMARCA4, SMARCE1)
- Histone modifications (H3K27me3, H3K4me3)
- SMARCB1 variants: WT, I315I (silent), I315R, I315X, W281P, W281X

**Dataset:**
- 72 samples total
- 6 genotypes × 5 targets × 2 replicates + IgG controls
- Sequenced on NovaSeq 6000
- ~30-50 million paired-end reads per sample

**Key Findings Enabled by Pipeline:**
- Nonsense mutants (I315X, W281X) show minimal SMARCB1 binding
- I315I (silent) retains WT-like binding patterns
- Missense mutants show partial loss of function
- SWI/SNF consensus regions validated therapeutic targets

## Troubleshooting Development Issues

### Pipeline Runs But No Output
**Cause:** Missing `publishDir` directive  
**Fix:** Add `publishDir "${params.outdir}/module_name", mode: 'copy'`

### Channels Empty / Processes Don't Execute
**Cause:** Filter removed all items or join failed  
**Fix:** Add `.view()` operators to debug channel contents

### Process Fails Silently
**Cause:** Error in bash script not caught by Nextflow  
**Fix:** Add `set -euo pipefail` at top of bash scripts, use `"""` for multi-line

### Parameters Not Applied
**Cause:** Parameter collision or wrong scope  
**Fix:** Use unique param names, check if workflow inherits correctly

### Resume Doesn't Work
**Cause:** Process signature changed (name, inputs, outputs)  
**Fix:** Don't change process names during development, or accept full re-run

### Memory Errors
**Cause:** Process memory config too low  
**Fix:** Increase in config: `withName: 'processName' { memory = '16 GB' }`

## Performance Optimization Notes

### What We Tried

**Symlink vs Copy:**
- Symlink: 6.5 min/sample
- Copy: 30+ min/sample
- **Winner:** Symlink (5× faster)

**Parallel Processing:**
- Sequential: ~8 hours for 72 samples
- Parallel (maxForks=18): ~5 hours
- **Winner:** Parallel with optimal fork count

**Channel Operations:**
- Multiple `.map().filter()` chains vs single complex operation
- Separate operations slightly slower but much easier to debug
- **Winner:** Readability over marginal performance gains

**BigWig Generation:**
- bamCoverage with multiple cores: 1-2 min/sample
- Splitting by chromosome: Overhead not worth it
- **Winner:** Single bamCoverage call with cpus=4

### Optimization Recommendations

1. **Use symlinks** for large files (already configured)
2. **Tune maxForks** based on available CPUs
3. **Profile with `-with-trace`** to identify bottlenecks
4. **Use local executor** for single machine (already configured)
5. **Don't over-parallelize** small/fast processes

## Project-Specific Notes

### Genotype Naming
- WT: Wild-type SMARCB1
- I315I: Silent variant (no protein change)
- I315R: Missense (conservative)
- I315X: Nonsense (no functional protein)
- W281P: Missense (disruptive)
- W281X: Nonsense (no functional protein)

### Antibodies/Targets
- SMARCB1, SMARCA4, SMARCE1: SWI/SNF core subunits
- H3K27me3: Polycomb repressive mark
- H3K4me3: Active promoter mark
- IgG: Negative control (background)

### Experimental Design
- 2 biological replicates per condition
- Paired-end sequencing (2×150 bp)
- Spike-in normalization (E. coli DNA)
- IgG controls for each genotype/replicate

### Data Location
- **AWS:** `/home/ec2-user/cutnrun-bam-processing/`
- **Processed bedgraphs:** `results/normalized_bedgraphs/`
- **Final peaks:** `results/filtered_peaks/`
- **Consensus regions:** `results/complex_binding/`
