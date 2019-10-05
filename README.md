
<!-- README.md is generated from README.Rmd. Please edit that file -->

# rgenesconverged

<!-- badges: start -->

<!-- badges: end -->

The goal of rgenesconverged is to identify instances of convergent
evolution by sequence similarity at the molecular level, by finding
conserved, similar, and independently arising sequences within
user-submitted phylogenetic trees. It’s modeled from concepts introduced
in a seminal paper:

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

## Example

This is a basic example which shows you how to solve a common problem:

``` r
#library(rgenesconverged)
# We use the example datasets from phangorn, tree and primates (see ./data).
load("~/rgenesconverged/data/smallTree.rda")
load("~/rgenesconverged/data/primates.rda")

summary(smallTree)
#>             Length Class  Mode     
#> edge        4      -none- numeric  
#> Nnode       1      -none- numeric  
#> tip.label   2      -none- character
#> edge.length 2      -none- numeric
summary(primates)
#>            Length Class  Mode   
#> Mouse      217    -none- numeric
#> Bovine     217    -none- numeric
#> Lemur      217    -none- numeric
#> Tarsier    217    -none- numeric
#> Squir Monk 217    -none- numeric
#> Jpn Macaq  217    -none- numeric
#> Rhesus Mac 217    -none- numeric
#> Crab-E.Mac 217    -none- numeric
#> BarbMacaq  217    -none- numeric
#> Gibbon     217    -none- numeric
#> Orang      217    -none- numeric
#> Gorilla    217    -none- numeric
#> Chimp      217    -none- numeric
#> Human      217    -none- numeric

#getConvergent(smallTree, primates, "Mouse", 1)
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />
