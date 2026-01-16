# Project Context & Development History

## Overview
This pipeline was developed over [timeframe] to automate CUT&RUN analysis from raw BAMs through peak calling and consensus analysis. Built for dissertation work on SMARCB1 variants in pediatric cancer, then generalized for broader use.

## Key Design Decisions

### Pipeline vs Analysis Separation
**Decision:** Main pipeline stops at replicate overlaps. Consensus analysis is a separate workflow.
**Rationale:** Clean separation between data processing (generalizable) and biological questions (project-specific). Keeps pipeline reusable while allowing custom downstream analysis.

### Adaptive Peak Filtering Strategy
**Decision:** Hybrid threshold = max(SNR from IgG, fixed threshold of 5)
**Rationale:** 
- Pure SNR fails with very clean IgG controls (everything passes)
- Pure fixed threshold ignores data quality
- Hybrid approach validated on nonsense mutants (I315X, W281X) that had no functional protein
- 58/60 samples used fixed threshold, 2 used higher SNR (I315R with noisy IgG)

### Why Not Include Master Lists in Main Pipeline?
**Decision:** Removed `createMasterPeakLists` that filtered for only WT/I315I
**Rationale:** That's analysis (selecting specific genotypes), not processing. Pipeline should create overlaps for ALL genotypes, then users decide what to compare.

### Nextflow Configuration Choices
- `maxForks = 18`: Parallelization sweet spot for 64 CPUs
- `stageInMode = 'symlink'`: Critical for BAM processing speed (5-10× faster)
- Separate param names (`consensus_outdir` vs `outdir`): Prevents conflicts when running multiple workflows

## Technical Implementation Notes

### Channel Operations Patterns
- Use `.map` to transform tuples
- Use `.filter` for subsetting
- Use `.join` to match on keys (with `by:` parameter)
- Use `.groupTuple` to collect by key
- Use `.combine` for Cartesian joins

### Parameter Handling
- Command-line comma-separated params need parsing: `params.genotypes.tokenize(',')`
- Always check `instanceof String` before tokenizing
- Validate parameters before workflow runs

### File Path Management
- `${projectDir}` for files in repo
- Absolute paths for external data
- Relative paths resolve from where `nextflow run` is executed
- Use unique param names to avoid conflicts

## Development Evolution

### Phase 1: Initial BAM Processing
Started with manual scripts, converted to modular Nextflow pipeline.

### Phase 2: Added Normalization
Spike-in normalization and bedgraph averaging.

### Phase 3: BigWig Generation
deepTools integration for visualization files.

### Phase 4: Peak Calling (Major Development)
**Initial approach:** Just run SEACR with fixed threshold
**Problem:** Nonsense mutants (no functional protein) still called thousands of peaks
**Solution 1:** Signal-to-noise ratio from IgG
**Problem with SNR:** Clean IgG controls = very low thresholds = everything passes
**Final solution:** Hybrid approach (max of SNR and fixed threshold)
**Validation:** I315X_SMARCB1 filtered from 109k peaks → 4 peaks (0.0% retention)

### Phase 5: Consensus Analysis
**Initial plan:** Include in main pipeline
**Pivoted:** Separate workflow for clean architecture
**Key insight:** SWI/SNF analysis is ONE use case, not THE use case

## Project-Specific Biological Context

### SMARCB1 Variants Studied
- **WT**: Wild-type (functional)
- **I315I**: Silent variant (functional, control)
- **I315R**: Missense, conservative substitution
- **I315X**: Nonsense, premature stop (NO functional protein)
- **W281P**: Missense, disruptive
- **W281X**: Nonsense, premature stop (NO functional protein)

### Why This Matters for Pipeline
Nonsense mutants provided critical validation of filtering strategy—they SHOULD have minimal peaks since there's no functional protein to bind chromatin.

### Dataset Characteristics
- 72 samples: 6 genotypes × 5 targets × 2 replicates + IgG controls
- Targets: SMARCB1, SMARCA4, SMARCE1 (SWI/SNF subunits), H3K27me3, H3K4me3 (histone marks)
- Platform: NovaSeq 6000, paired-end
- Processing: AWS EC2 m7i.16xlarge (64 CPUs)

## Key Learnings

### What Worked Well
1. Modular architecture = easy to iterate on individual components
2. `-resume` flag = saved hours during development
3. Separate analysis workflows = keeps pipeline clean
4. Hybrid filtering = validated on real biological negatives

### What Didn't Work
1. Pure SNR filtering = too lenient with clean IgG
2. Including project-specific analysis in main pipeline = poor architecture
3. Not using symlinks initially = 5× slower

### Development Gotchas
- **bedtools warnings**: "Query marked as paired but mate does not occur next to it" = SAFE TO IGNORE
- **Process config warnings**: When running consensus workflow, warnings about missing processes = SAFE TO IGNORE
- **Resume breaks**: Changing process names breaks cache
- **publishDir**: Always use unique param names to avoid conflicts
- **Parameter inheritance**: Workflows inherit params from main config

## Future Development Plans

### Immediate Next Steps
1. Automated spike-in alignment (BWA) and normalization factor calculation
2. Data visualization module (heatmaps, profiles)
3. Quality control report (MultiQC integration)

### Medium Term
1. Differential binding analysis (DiffBind)
2. Peak annotation (ChIPseeker)
3. Expand to ATAC-seq, RNA-seq, variant calling

### Long Term
1. Containerization (Docker/Singularity)
2. nf-core style modularization
3. Web interface for non-command-line users

## Important Context for Future Development

### When to Use Main Pipeline
- Processing new CUT&RUN data from BAMs
- Need standardized, reproducible peak calls
- Want QC metrics and replicate concordance

### When to Create New Analysis Workflows
- Have specific biological question
- Need custom peak comparisons
- Want different genotype/target combinations

### Architecture Principles to Maintain
1. **Main pipeline = data processing** (BAMs → peaks)
2. **Analysis workflows = biological questions** (peak comparisons)
3. **Keep modules independent** (can run separately if needed)
4. **Document everything** (future you will thank you)

## Performance Benchmarks (for reference)

| Operation | Samples | Time | Notes |
|-----------|---------|------|-------|
| BAM processing | 1 | 6.5 min | Bottleneck |
| Peak calling | 60 | 45 min | Parallel |
| Full pipeline | 72 | ~5 hours | End-to-end |

## Testing Strategy

### Quick Test
```bash
nextflow run main.nf --sample_list samples.txt --test_mode true
```

### Full Validation
```bash
# Run on known good samples
# Check nonsense mutants get filtered
# Verify replicate concordance
# Compare to manual analysis results
```

## Contact & Maintenance

**Developer:** Garrett Cooper (garrettc00per)
**Institution:** Emory University
**Status:** PhD candidate → Industry transition
**Last Major Update:** [Current date]

## References to Key Discussions

### Filtering Strategy Evolution
See commit history for discussion of SNR vs fixed vs hybrid approaches.

### Pipeline Architecture Decisions
Main pipeline vs analysis workflow separation discussed extensively—consensus was to keep processing separate from biological analysis.

### Parameter Naming Conventions
Decision to use unique names (`consensus_outdir`) to avoid conflicts documented in DEVELOPMENT.md.

## Questions This Document Should Answer

1. Why hybrid filtering? → Validated on nonsense mutants
2. Why separate consensus workflow? → Architecture principle
3. Why these specific parameters? → Optimized through testing
4. What's next? → See future development plans
5. How to extend? → Follow modular architecture pattern

---

*This document should be updated whenever major design decisions are made or significant features are added.*
