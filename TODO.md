# TODO

## Completed ✅
- [x] BAM Processing Module
- [x] Bedgraph Normalization Module
- [x] BigWig Generation Module
- [x] Test full 72-sample run

## Immediate Next Steps

### Peak Calling Module (Priority 1)
- [ ] Create modules/peak_calling.nf
- [ ] Implement SEACR for CUT&RUN data
- [ ] Add stringent and relaxed thresholds
- [ ] Use normalized bedgraphs + IgG controls
- [ ] Filter peaks by minimum intensity (>5)
- [ ] Generate reproducible peaks (R1 ∩ R2 overlap)

### Spike-in Alignment Module (Priority 2)
- [ ] Create modules/spike_in_normalization.nf
- [ ] Align FASTQs to human (GRCh38) + E.coli
- [ ] Remove duplicates with Picard
- [ ] Calculate normalization factors automatically
- [ ] Output: normalization_factors.tsv

### Visualization & QC (Priority 3)
- [ ] Add bedGraph IgG subtraction for averaged files
- [ ] Create scatter plot generation (Mutant vs WT)
- [ ] Generate correlation heatmaps
- [ ] Signal distribution plots

## Future Modules

### Differential Binding Analysis
- [ ] DiffBind or csaw integration
- [ ] Condition comparisons
- [ ] Statistical testing
- [ ] Volcano plots and MA plots

### Documentation
- [ ] Add peak calling usage examples
- [ ] Document BigWig visualization workflows
- [ ] Create troubleshooting guide for deepTools

## Questions to Resolve
- Peak calling: SEACR stringent vs relaxed threshold?
- Should we filter averaged bedGraphs by peak regions only?
- Add IGV session file generation?
