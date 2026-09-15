# Transform a list of Pathway IDs into a binary matrix.

This function transforms a list of \\p\\ (one vector of pathway IDs per
gene) groups into a binary matrix.

## Usage

``` r
target.help(genes)
```

## Arguments

- genes:

  List of \\p\\ items. Each item is the vector of Pathway IDs a gene
  belongs to.

## Value

A \\p \times p\\ binary matrix: the coefficient (i,j) is 1 if genes i
and j belong to a common pathway and 0 otherwise.

## See also

[`targetF`](https://github.com/vguillemot/SHIP/reference/targetF.md),[`targetG`](https://github.com/vguillemot/SHIP/reference/targetG.md),[`targetGpos`](https://github.com/vguillemot/SHIP/reference/targetGpos.md),
[`targetGstar`](https://github.com/vguillemot/SHIP/reference/targetGstar.md).

## Author

Monika Jelizarow and Vincent Guillemot

## Examples

``` r

g1 <- c("path1", "path2", "path3", "path4")
g2 <- c("path5", "path6", "path3", "path11")
g3 <- c("path10", "path5", "path12", "path13")
target.help(list(g1, g2, g3)) 
#>      [,1] [,2] [,3]
#> [1,]    1    1    0
#> [2,]    1    1    1
#> [3,]    0    1    1

```
