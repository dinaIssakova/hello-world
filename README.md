
<!-- README.md is generated from README.Rmd. Please edit that file -->

# rgenesconverged

<!-- badges: start -->

<!-- badges: end -->

The goal of rgenesconverged is to identify instances of convergent
evolution by sequence similarity at the molecular level, by finding
conserved, similar, and independently arising sequences within
user-submitted phylogenetic trees. For the mathematical model on which
this package is based, please see:

Zhang, J. and Kumar, S. (1997) Detection of Convergent and Parallel
Evolution at the Amino Acid Sequence Level. Mol. Biol. Evol.
14(5):527-536.

## Installation

You can install the first-draft version of rgenesconverged from github
with:

``` r
library(devtools)
install_github("dinaIssakova/rgenesconverged")
library(rgenesconverged)
```

## Overview

An overview of the package is illustrated below.

![](./inst/extdata/ISSAKOVA_D_A1.png)

## Contributions

The author of the package is Dina Issakova. The functions available
within this package include:

    rlibrary("rgenesconverged")
    lsf.str("package:rgenesconverged")

  - getBranchLength
  - getMostRecentCommonAncestor
  - mapLetters
  - areCondSatisfied
  - probOfChange
  - probOfSiteConfig
  - probOfNSitesByChance
  - getConvergent
  - convSiteData
  - rgenesconvergedPlot

All functions above were authored by Dina Issakova. Examples and help
files make use of datasets available from phangorn R package and subsets
of these datasets. PAM matrix included in dataset is manually copied
from:

Dayhoff, M.O. et al. (1978) A model of evolutionary change in proteins.
Atlas of Protein Sequence and Structure. p. 345-352

The getMostRecentCommonAncestor function makes use of commonly known
algorithms for finding the deepest ancestor of two leaves in a binary
tree. The areCondSatisfied and probOfSiteConfig functions make use of
the phangorn R package function ancestral.pars to predict ancestral
states of amino acid sequences based on branch length and mutation
probabilities. The areCondSatisfied, probOfChange and probOfSiteConfig
functions make use of the Biostrings R package to translate amino acid
sequences using the standard genetic code and access standard lists of
amino acids.

The areCondSatisfied, probOfChange, probOfSiteConfig,
probOfNSitesByChance, getConvergent and convSiteData implement
mathematical models of genomic convergent evolution described in :

Zhang, J. and Kumar, S. (1997) Detection of Convergent and Parallel
Evolution at the Amino Acid Sequence Level. Mol. Biol. Evol.
14(5):527-536.

The rgenesconvergedPlot was authored by Dina Issakova and makes use of
the package ggtree to plot the phylogenetic tree.
