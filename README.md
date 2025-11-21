# Genome-Host Association (GHA) Analysis
## Genome-host association approach for fungal pathogen genes involved in host specialisation
This guide outlines the full GHA workflow from raw sequencing reads to candidate gene identification and downstream functional validation. It is designed for researchers investigating host specialization or cultivar-specific virulence factors in plant pathogens.

This repository gathers the scripts that I have used to explore the genetic basis of host specialization in the wheat fungal pathogen *Zymoseptoria tritici*. The biorXiv version of the associated paper will be added. The link will be updated upon acceptance of the publication. 

## Step 1: Data preparation

SNP calling and filtering were performed as described in Feurtey et al., 2023 [[WW_PopGen](https://github.com/afeurtey/WW_PopGen?tab=readme-ov-file#data-preparation)].

## Step 2: Genome-host association mapping (GHA)
We performed a genome-host association (GHA) study based on GWAS using the host cultivar as a phenotype with a “one-versus-all” approach. We implemented a categorical GWAS by comparing strains sampled from a given wheat cultivar (phenotype 1) versus strains from the other 11 cultivars (phenotype 0). To avoid unbalanced sample sizes, each given group of samples (e.g. 64 strains from the cultivar Arina) was compared to an equal number of samples randomly selected from the rest of the collection. We repeated this random sub-sampling step 100 times. The phenotype matrix per cultivar can be found [[here](Matrix_host_phenotypes)].

We applied this approach across all 12 hosts, generating a total of 1200 GWAS analyses run using [[vcf2gwas](https://github.com/frankvogt/vcf2gwas)], with a univariate linear mixed model Wald test. The script for GHA is available [[here](GHA)], with example files. Significant SNPs were corrected with the Bonferroni method by dividing 0.05 by the sum of unique variants analyzed.

The script for GHA visualization can be found [[here](GHA)].

## Functional validation of GHA candidate genes
The plasmids maps used for *Z. tritici* transformations in Genebank format are available [[here](plasmids_maps)].
