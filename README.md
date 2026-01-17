# NGS Project - Group 15 - PRJNA340216 (Diet vs gut microbiome)

**Dataset**: [PRJNA340216](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA340216)

**Main reference**: De Angelis, M., Ferrocino, I., Calabrese, F.M. et al. Diet influences the functions of the human intestinal microbiome. Sci Rep 10, 4247 (2020). https://doi.org/10.1038/s41598-020-61192-y

**Scope**: Start with 2 diet groups (vegans vs omnivores) and benchmark runtime on a small subset.

## Project layout

```
- docs/      notes + plan
- scripts/   helper scripts
- workflow/  pipeline definitions (optional)
- data/      NOT tracked (raw/ processed/ reference/)
- results/   figures/tables (small outputs only)
- logs/      small logs
```

## Running scripts

### Preprocessing

To preprocess single sample, run:

```bash
sbatch /home/projects/22126_NGS/projects/group15/scripts/hpc/run_preprocess.sh /home/projects/22126_NGS/projects/group15/data/raw/mysample.fastq.gz
```

If a target output file exits, the corresponding preprocessing step (QC, trimming, decontamination) will be skipped. **Note**: it's a very crude check for presence of expected output file by name, it does not evaluate whether the input file content has changed or similar.
