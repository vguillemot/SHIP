
<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- badges: start -->

<a href="https://cran.r-project.org/package=SHIP"
class="SHIP-release"><img
src="https://www.r-pkg.org/badges/version/SHIP" alt="CRAN Status" /></a>
<!-- badges: end -->

# SHrinkage covariance Incorporating Prior knowledge

<!-- badges: start -->

<!-- badges: end -->

The SHIP package implements the shrinkage estimator of a covariance
matrix given any covariance target, such as described by Schaefer and
Strimmer in 2005. In addition, it proposes several targets based on
biological knowledge extracted from the public database KEGG.

To use the shrinkage estimator, one should just have at hand a data set
in the form of a $n \times p$ matrix, and a covariance target.

If one wishes to use the proposed targets, the data set should be
compatible with KEGG, i.e. it should be possible to extract for each
gene the pathways it belongs to. This information, for example, can be
found in libraries such as hgu133plus2.db. \## Installation

You can install the development version of SHIP from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("vguillemot/SHIP")
```
