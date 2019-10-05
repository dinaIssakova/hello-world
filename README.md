
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
## basic example code
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
