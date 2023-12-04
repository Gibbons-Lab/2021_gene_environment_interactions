# Material and code for the Genome Interplay Manuscript

This repository contains materials and code for the manuscript:

*Genome-microbiome interplay provides insight into the determinants of the human blood metabolome*<br>
Christian Diener#, Chengzhen L. Dai#, Tomasz Wilmanski, Priyanka Baloni, Brett Smith, Noa Rappaport, Leroy Hood, Andrew T. Magis and Sean M. Gibbons<br>
https://doi.org/10.1101/2022.02.04.479172

## Organization of the repository

```
root
    > notebook.ipynb  # Analysis steps in Jupyter notebooks
    > workflow.nf     # Nextflow workflow for larger analyses
    > figures
        + fig1.png    # Figure 1 of the manuscript
        + fig1.svg    # SVG source file for multi-panel figures
        [...]
    > data
        + data.csv    # Intermediate data files or results tables
        [...]
    [...]
```

## Installation of dependencies

All dependencies can be installed with the provided [conda](https://docs.conda.io/en/latest/miniconda.html) environment.

```bash
conda env create -f conda.yml
```

This will create an environment `gxe2021` that can be used to run all analyses.

```bash
conda activate gxe2021
```

## Analysis steps

The analyses from the manunscript can be reproduced with the following steps. This will require obtaining the data from the Arivale cohort first.

> Qualified researchers can access the full Arivale deidentified dataset supporting the findings in this study for research purposes through signing a Data Use Agreement (DUA). Inquiries to access the data can be made at data-access{at}isbscience.org and will be responded to within 7 business days.

Raw sequencing data for the 16S amplicon sequencing has been deposited to the SRA. A full list of run accessions can be found in the `data` directory [here](https://github.com/Gibbons-Lab/2021_gene_environment_interactions/blob/main/data/sra_manifest.csv).

### Analysis steps

1. Assembly of the cohort | [notebook](cohort.ipynb)
2. Confounder adjustment  | [notebook](residuals.ipynb)
3. Microbiome-metabolite associations for the training cohort | [notebook](microbe_associations.ipynb)
4. mGWAS for the training cohort | [workflow](gwas.nf)<br>
   Run this with `nextflow run -resume gwas.nf`
5. Inspect results from the mGWAS | [notebook](gwas.ipynb)
6. Fit models and obtain R2 for training cohort | [notebook](train_variances.ipynb)
7. Obtain out-of-sample R2 for validation cohort | [notebook](validation_variances.ipynb)
8. Detailed analysis of bile acid and sphingolipid R2 | [notebook](BAs_and_sphingolipids.ipynb)
9. Fit genome-microbiome-metabolome interactions | [notebook](genome_microbiome_interactions.ipynb)
10. Analyze the obtained interactions | [notebook](inspect_genome_microbiome_interactions.ipynb)

