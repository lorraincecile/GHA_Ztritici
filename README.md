# Genome-Host Association (GHA) Analysis
## Genome-host association approach for fungal pathogen genes involved in host specialisation
This guide outlines the full GHA workflow from raw sequencing reads to candidate gene identification and downstream functional validation. It is designed for researchers investigating host specialization or cultivar-specific virulence factors in plant pathogens.

This repository gathers the scripts that I have used to explore the genetic basis of host specialization in the wheat fungal pathogen *Zymoseptoria tritici*. The biorXiv version of the associated paper will be added. The link will be updated upon acceptance of the publication. 

## Step 1: Data preparation

SNP calling and filtering were performed as described in Feurtey et al., 2023 [[WW_PopGen](https://github.com/afeurtey/WW_PopGen?tab=readme-ov-file#data-preparation)].

## Step 2: Genome-host association mapping (GHA)

We performed a GHA study using a “one-versus-all” categorical GWAS framework.  
For each host cultivar:

- `1` = strains collected from the focal host  
- `0` = strains collected from any of the other 11 hosts  

To avoid unbalanced sample sizes and maintain statistical power, each host group (e.g., 64 strains from cultivar Arina) was compared to the **same number of randomly sampled strains** drawn from the remainder of the population.  
This **random subsampling was performed 100 times per host**.

The resulting phenotype matrices are available here:  
[Matrix_host_phenotypes](Matrix_host_phenotypes)

### GHA analyses using `vcf2gwas`

Across 12 cultivars × 100 subsamples, a total of **1200 GWAS scans** were performed using  
[`vcf2gwas`](https://github.com/frankvogt/vcf2gwas) with a univariate LMM (GEMMA).

#### Script used to run GHA

```bash
vcf=data/vcf/filtered_dataset.vcf.gz
input=data/phenotypes/<host_name>/

for subsample_file in ${input}*.csv
do
    phenotypes=$(basename ${subsample_file})
    output=$(basename ${subsample_file%.csv})

    echo "Running GHA for: ${phenotypes}"

    vcf2gwas -v ${vcf} \
             -pf ${input}${phenotypes} \
             -ap \
             -lmm \
             -nl \
             -np \
             -o ${input}${output}/
done
```
### Negative-control GHA tests (mock and null)

To empirically evaluate the tendency of the GHA workflow to produce false-positive host associations, we ran two negative-control analyses where group labels are randomized while keeping the genotype data, SNP filtering, covariates, and LMM-GWAS settings identical to the main analyses.

- **Mock GHA (pseudo focal group vs re-sampled control group):** randomly select 80 isolates as a pseudo “focal host” group (independent of host-of-origin) and compare them to the remaining isolates. Repeat the full bootstrap GWAS procedure (100 iterations; resample the 80 isolates each time).  
  **Outcome:** no SNP reached the Bonferroni threshold.

- **Null GHA (random vs random):** randomly sample two groups of equal size (80 vs 80 isolates), independent of host-of-origin, and run GWAS with the same pipeline.  
  **Outcome:** no SNP reached the Bonferroni threshold.
  

### Visualization GHA
The script for GHA visualization can be found [[here](GHA)].

## Functional validation of GHA candidate genes
In this study, we prioritized **functional validation** of GHA candidates via a targeted deletion screen in *Z. tritici* as an empirical “positive control” for candidate relevance. More generally, we recommend that future applications of this pipeline also include at least one **a priori positive-control association** (e.g., a cultivar/host factor with a known genotype expected to map to a specific pathogen locus) alongside mock/null negative controls to benchmark sensitivity and interpretability.

The plasmids maps used for *Z. tritici* transformations in Genebank format are available [[here](plasmids_maps)].
