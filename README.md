# RVS - Rare Variant Sharing R Package

[![Bioc](https://bioconductor.org/images/logo_bioconductor.gif)](https://bioconductor.org/packages/RVS)
[![downloads](https://bioconductor.org/shields/downloads/release/RVS.svg)](http://bioconductor.org/packages/stats/bioc/RVS/)
[![Build Status](https://travis-ci.org/sherman5/RVS.svg?branch=master)](https://travis-ci.org/sherman5/RVS)

Rare Variant Sharing (RVS) implements tests of association and linkage between rare genetic variant genotypes and a dichotomous phenotype, e.g. a disease status, in family samples. The tests are based on probabilities of rare variant sharing by relatives under the null hypothesis of absence of linkage and association between the rare variants and the phenotype and apply to single variants or multiple variants in a region (e.g. gene-based test).

# Installing RVS

*RVS* is a bioconductor R package and so the release version can be installed
as follows:

```
install.packages("BiocManager")
BiocManager::install("RVS")
```

The most up-to-date version of *RVS* can be installed directly from Github:

```
## Method 1 using BiocManager
BiocManager::install("sherman5/RVS")

## Method 2 using devtools package
devtools::install_github("sherman5/RVS")
```

# Using RVS

Follow the vignette here: https://bioconductor.org/packages/release/bioc/vignettes/RVS/inst/doc/RVS.html

# Reporting Bugs/Getting Help

If you encounter a bug in the package, or have any questions about RVS, please open an issue at https://github.com/sherman5/RVS/issues
