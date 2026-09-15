# A short introduction to SHIP

SHIP estimates covariance matrices by combining the sample covariance
matrix with a structured target. The final estimate has the form

\\ \widehat{\Sigma} = \lambda T + (1 - \lambda)S, \\

where \\S\\ is the sample covariance matrix, \\T\\ is a target, and \\0
\leq \lambda \leq 1\\ controls the amount of shrinkage.

## The example data

The package includes `expl`, a small microarray data set with 102
observations and 100 genes. The accompanying `genegroups` list contains
KEGG pathway IDs when they are available; `NA` means that no pathway
information is available for that gene.

Code

\
` `[`data`](https://rdrr.io/r/utils/data.html)`(``"expl"``)`\
[`dim`](https://rdrr.io/r/base/dim.html)`(``expl``$``x``)`\
`#> [1] 102 100`\
[`length`](https://rdrr.io/r/base/length.html)`(``expl``$``genegroups``)`\
`#> [1] 100`\
[`sum`](https://rdrr.io/r/base/sum.html)`(`[`vapply`](https://rdrr.io/r/base/lapply.html)`(``expl``$``genegroups``, ``function``(``group``)`` ``!`[`all`](https://rdrr.io/r/base/all.html)`(`[`is.na`](https://rdrr.io/r/base/NA.html)`(``group``)``)``, `[`logical`](https://rdrr.io/r/base/logical.html)`(``1``)``)``)`\
`#> [1] 37`

## Comparing covariance targets

[`build.target()`](https://github.com/vguillemot/SHIP/reference/build.target.md)
constructs several targets. `D` and `F` use only the data; `G`, `Gpos`,
`Gstar`, and `cor` also use the pathway groups.

Code

\
`target_types`` ``<-`` `[`c`](https://rdrr.io/r/base/c.html)`(``"D"``, ``"F"``, ``"G"``, ``"Gpos"``, ``"Gstar"``, ``"cor"``)`\
`targets`` ``<-`` `[`setNames`](https://rdrr.io/r/stats/setNames.html)`(`\
`  `[`lapply`](https://rdrr.io/r/base/lapply.html)`(``target_types``, ``function``(``type``)`` ``{`\
`    `[`build.target`](https://github.com/vguillemot/SHIP/reference/build.target.md)`(``expl``$``x``, ``expl``$``genegroups``, type ``=`` ``type``)`\
`  ``}``)``,`\
`  ``target_types`\
`)`\
\
`target_summary`` ``<-`` `[`data.frame`](https://rdrr.io/r/base/data.frame.html)`(`\
`  target ``=`` `[`names`](https://rdrr.io/r/base/names.html)`(``targets``)``,`\
`  dimension ``=`` `[`vapply`](https://rdrr.io/r/base/lapply.html)`(``targets``, ``function``(``target``)`` `[`paste`](https://rdrr.io/r/base/paste.html)`(`[`dim`](https://rdrr.io/r/base/dim.html)`(``target``)``, collapse ``=`` ``" x "``)``, `[`character`](https://rdrr.io/r/base/character.html)`(``1``)``)``,`\
`  nonzero ``=`` `[`vapply`](https://rdrr.io/r/base/lapply.html)`(``targets``, ``function``(``target``)`` `[`sum`](https://rdrr.io/r/base/sum.html)`(``target`` ``!=`` ``0``)``, `[`integer`](https://rdrr.io/r/base/integer.html)`(``1``)``)``,`\
`  row.names ``=`` ``NULL`\
`)`\
`target_summary`\
`#>   target dimension nonzero`\
`#> 1      D 100 x 100     100`\
`#> 2      F 100 x 100   10000`\
`#> 3      G 100 x 100    4028`\
`#> 4   Gpos 100 x 100    2252`\
`#> 5  Gstar 100 x 100    4028`\
`#> 6    cor 100 x 100    3012`

The diagonal target keeps the sample variances and sets covariances to
zero. The constant-correlation target uses one average correlation for
all off-diagonal entries. The knowledge-based targets restrict
correlations to genes that share a pathway; `Gpos` keeps positive links,
while `Gstar` separates positive and negative links.

## Estimating a covariance matrix

The target can then be passed to
[`shrink.estim()`](https://github.com/vguillemot/SHIP/reference/shrink.estim.md).

Code

\
`estimate`` ``<-`` `[`shrink.estim`](https://github.com/vguillemot/SHIP/reference/shrink.estim.md)`(``expl``$``x``, ``targets``[[``"G"``]``]``)`\
\
[`c`](https://rdrr.io/r/base/c.html)`(`\
`  lambda ``=`` ``estimate``$``lambda``,`\
`  sample_variance ``=`` `[`sum`](https://rdrr.io/r/base/sum.html)`(`[`diag`](https://rdrr.io/r/base/diag.html)`(`[`cov`](https://rdrr.io/r/stats/cor.html)`(``expl``$``x``)``)``)``,`\
`  estimated_variance ``=`` `[`sum`](https://rdrr.io/r/base/sum.html)`(`[`diag`](https://rdrr.io/r/base/diag.html)`(``estimate``$``shrink.cov``)``)`\
`)`\
`#>             lambda    sample_variance estimated_variance `\
`#>         0.04397249        24.93098682        24.93098682`

A value of `lambda` close to zero keeps more of the empirical covariance
matrix. A value close to one gives more weight to the structured target.

The target definitions are based on Schaefer and Strimmer (2005) and
Jelizarow et al. (2010), *Over-optimism in bioinformatics: an
illustration*.
