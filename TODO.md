# TODO

## Immediate Next Steps
- [ ] Test full 72-sample run
- [ ] Verify all output files are correct
- [ ] Add stageInMode to nextflow.config if not present

## Future Modules

### Spike-in QC Module (Priority 1)
- [ ] Create modules/spike_in_qc.nf
- [ ] Calculate % aligned to E.coli vs human
- [ ] Output CSV with stats per sample
- [ ] Create bin/calculate_normalization.py script

### Peak Calling Module (Priority 2)
- [ ] Research SEACR parameters for CUT&RUN
- [ ] Create modules/peak_calling.nf
- [ ] Integrate optional normalization
- [ ] Test on subset of samples

### Documentation
- [ ] Write detailed spike-in normalization guide
- [ ] Add usage examples for each module
- [ ] Create troubleshooting section

## Questions to Resolve
- What SEACR threshold to use?
- How to handle IgG controls in peak calling?
- Should differential binding be included or separate pipeline?
