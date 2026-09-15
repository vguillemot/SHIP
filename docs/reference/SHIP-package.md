# SHIP provides shrinkage covariance estimation with user-selected targets. The available targets include diagonal, constant-correlation, and knowledge-based structures informed by functional gene groups.

Start with
[`build.target`](https://github.com/vguillemot/SHIP/reference/build.target.md)
to construct a target matrix, then pass it to
[`shrink.estim`](https://github.com/vguillemot/SHIP/reference/shrink.estim.md)
together with the data matrix.

## References

- J. Schaefer and K. Strimmer, 2005. A shrinkage approach to large-scale
  covariance matrix estimation and implications for functional genomics.
  Statist. Appl. Genet. Mol. Biol. 4:32.

- M. Jelizarow, V. Guillemot, A. Tenenhaus, K. Strimmer, A.-L.
  Boulesteix, 2010. Over-optimism in bioinformatics: an illustration.
  Bioinformatics. Accepted.

## See also

Useful links:

- <https://github.com/vguillemot/SHIP>

- Report bugs at <https://github.com/vguillemot/SHIP/issues>

## Author

Monika Jelizarow and Vincent Guillemot

## Examples

``` r
data("expl")
target <- build.target(expl$x, expl$genegroups, type = "G")
estimate <- shrink.estim(expl$x, target)
estimate$lambda
#> [1] 0.04397249
```
