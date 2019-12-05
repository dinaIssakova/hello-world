
<!-- README.md is generated from README.Rmd. Please edit that file -->

# rgenesconverged

<!-- badges: start -->

<!-- badges: end -->

Sequence convergence is a type of convergent evolution that results in
similarity between orthologous genetic sequences. However, this can be
difficult to evaluate statistically given classical paradigms for
convergent evolution, which rest on the assumption that convergent
phenotypes are achieved by unrelated genetic mechanisms. While
sophisticated models exist to model this process and to evaluate the
probability of this explanation of genetic similarity across two
different species, no tool has currently been implemented in R to carry
out these analyses. Here we present rgenesconverged, an R package
allowing the exploration of sequence convergence hypotheses in a given
phylogeny.

This package implements and extends the model for sequence convergence
given by Zhang and Kumar:

J. Zhang and S. Kumar, “Detection of convergent and parallel evolution
at the amino acidsequence level,” Mol. Biol. Evol., vol. 14, no. 5,
pp. 527–536, 1997.

For more information or to use rgenesconverged in a publication, please
consult/cite:

Issakova, D. 2019. rgenesconverged: An R Package for the Exploration of
Molecular Convergent Evolution. bioRxiv doi:
<https://doi.org/10.1101/858076>

## Installation

You can install the current version of rgenesconverged from github with:

``` r
library(devtools)
install_github("dinaIssakova/rgenesconverged")
library(rgenesconverged)
```

And run the Shiny app by:

``` r
shinyRgenesconverged()
```

## Overview

``` r
browseVignettes(rgenesconverged)
```

The inputs required to rgenesconverged are an R phydat object containing
alignment information for the selected species and an R tree object - in
the RShiny app, only the alignment is asked for and the tree is built
from that. The user selects a species of interest for comparison to
which the other species in the tree are compared. Sequences are aligned
by ClustalOmega to the reference sequence; the position at which
convergence is to be evaluated corresponds to that position in the
reference sequence, and corresponding positions in the alignment. Two
possible ‘paths’ exist in this package for phylogenetic analysis. The
first is a direct implementation of Zhang and Kumar’s model, which is a
well-established and commonly used model for sequence convergence; the
second is a natural extension allowing for a more nuanced and
exploratory analysis, but less statistical rigor. An R plot is also
available, which finds species that are convergent at the particular
position specified to the reference species, and highlights them in a
phylogram. While designed for use in an R environment, an RShiny UI is
available with the package for users unfamiliar with R programming.

For more information, please see the vignettes for this package, or

Issakova, D. 2019. rgenesconverged: An R Package for the Exploration of
Molecular Convergent Evolution. bioRxiv doi:
<https://doi.org/10.1101/858076>

![](./inst/extdata/ISSAKOVA_D_A1.png)

## Contributions

    rlibrary("rgenesconverged")
    lsf.str("package:rgenesconverged")

The author of the package is Dina Issakova. The functions available
within this package include:

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
  - shinyRgenesconverged

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
amino acids. The probOfChange function makes use of the package foreach
for parallelizing.

The areCondSatisfied, probOfChange, probOfSiteConfig,
probOfNSitesByChance, getConvergent and convSiteData implement
mathematical models of genomic convergent evolution described in :

Zhang, J. and Kumar, S. (1997) Detection of Convergent and Parallel
Evolution at the Amino Acid Sequence Level. Mol. Biol. Evol.
14(5):527-536.

The rgenesconvergedPlot was authored by Dina Issakova and makes use of
the package ggtree and ggplot to plot the phylogenetic tree.

Specifically, Dina Issakova’s contribution was 1) implementing the
mathematical model described by Zhang and Kumar above and 2) extending
the model by adding a score-based exploratory function.

## Walkthrough

For walkthrough, please see vignettes.
