# Diagonal covariance target.

The target keeps the sample variances and sets all covariances to zero.
If \\S = (s\_{ij})\\ is the sample covariance matrix, the target is
\$\$T\_{ij} = \begin{cases} s\_{ii} & \text{if } i = j \\ 0 &
\text{otherwise} \end{cases}\$\$

## Usage

``` r
targetD(x, genegroups = NULL)
```

## Arguments

- x:

  A \\n \times p\\ data matrix.

- genegroups:

  The genegroups are not used for this target.

## Value

A \\p \times p\\ diagonal matrix.

## See also

Other covariance targets:
[`targetCor()`](https://github.com/vguillemot/SHIP/reference/targetCor.md),
[`targetF()`](https://github.com/vguillemot/SHIP/reference/targetF.md),
[`targetG()`](https://github.com/vguillemot/SHIP/reference/targetG.md),
[`targetGpos()`](https://github.com/vguillemot/SHIP/reference/targetGpos.md),
[`targetGstar()`](https://github.com/vguillemot/SHIP/reference/targetGstar.md)

## Author

Monika Jelizarow and Vincent Guillemot

## Examples

``` r

x <- matrix(rnorm(10*30),10,30)
target_matrix <- targetD(x, NULL)
```
