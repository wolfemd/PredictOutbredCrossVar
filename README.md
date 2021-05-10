# Reproducible documentation: Wolfe et al. Genomic mating in outbred species: predicting cross usefulness with additive and total genetic covariance matrices

This repository and website documents all analyses, summary, tables and figures associated with the following PREPRINT: [Genomic mating in outbred species: predicting cross usefulness with additive and total genetic covariance matrices](https://doi.org/10.1101/2021.01.05.425443)!

A [workflowR](https://workflowr.io/) project.

# PREPRINT

This repository and website documents all analyses, summary, tables and figures associated with the following PREPRINT: [Genomic mating in outbred species: predicting cross usefulness with additive and total genetic covariance matrices](https://doi.org/10.1101/2021.01.05.425443)!

# Abstract

Diverse crops are both outbred and clonally propagated. Breeders typically use truncation selection of parents and invest significant time, land and money evaluating the progeny of crosses to find exceptional genotypes. We developed and tested genomic mate selection criteria suitable for organisms of arbitrary homozygosity level where the full-sibling progeny are of direct interest as future parents and/or cultivars. We extended cross variance and covariance variance prediction to include dominance effects and predicted the multivariate selection index genetic variance of crosses based on haplotypes of proposed parents, marker effects and recombination frequencies. We combined the predicted mean and variance into usefulness criteria for parent and variety development. We present an empirical study of cassava (Manihot esculenta), a staple tropical root crop. We assessed the potential to predict the multivariate genetic distribution (means, variances and trait covariances) of 462 cassava families in terms of additive and total value using cross-validation. Most variance (89%) and covariance (70%) prediction accuracy estimates were greater than zero. The usefulness of crosses were accurately predicted with good correspondence between the predicted and the actual mean performance of family members breeders selected for advancement as new parents and candidate varieties.  We also used a directional dominance model to quantify significant inbreeding depression for most traits. We predicted 47,083 possible crosses of 306 parents and contrasted them to those previously tested to show how mate selection can reveal new potential within the germplasm. We enable breeders to consider the potential of crosses to produce future parents (progeny with top breeding values) and varieties (progeny with top own performance). 

# Start here

[**Project Homepage**](https://wolfemd.github.io/PredictOutbredCrossVar/)

# Data availability and reproducibility

The R package **workflowr** was used to document this study reproducibly.

Much of the supporting data *and* output from the analyses documented here are too large for GitHub.

The repository will be mirrored, here: <ftp://ftp.cassavabase.org/manuscripts/Wolfe_et_al_2021>

or until publication [here](ftp://ftp.cassavabase.org/marnin_datasets/).

# Supporting R package `predCrossVar`

In addition, we combined many of the core (and useful) support functions for predicting crosses into an R package **predCrossVar**, which is available on GitHub and can be installed with, e.g.:

```r
devtools::install_github("wolfemd/predCrossVar", ref = 'master') 
```

The functions in **predCrossVar** are used throughout.

# Supporting functions `code/`

The analyses in the **html** / **Rmd** files referenced above often source R scripts in the `code/` sub-folder. These are wrapper functions around the packaged core functions in **predCrossVar**, to do the specific analyses for this paper.
