# Package index

## Main functions

- [`build.target()`](https://github.com/vguillemot/SHIP/reference/build.target.md)
  : Creating a covariance target, optionally by using information from
  KEGG pathways.
- [`shrink.estim()`](https://github.com/vguillemot/SHIP/reference/shrink.estim.md)
  : Shrinkage estimator of the covariance matrix, given a data set and a
  covariance target.

## Covariance targets

- [`targetD()`](https://github.com/vguillemot/SHIP/reference/targetD.md)
  : Diagonal covariance target.
- [`targetF()`](https://github.com/vguillemot/SHIP/reference/targetF.md)
  : Constant-correlation covariance target.
- [`targetCor()`](https://github.com/vguillemot/SHIP/reference/targetCor.md)
  : Significance-filtered knowledge-based target.
- [`targetG()`](https://github.com/vguillemot/SHIP/reference/targetG.md)
  : Knowledge-based constant-correlation target.
- [`targetGpos()`](https://github.com/vguillemot/SHIP/reference/targetGpos.md)
  : Positive-correlation knowledge-based target.
- [`targetGstar()`](https://github.com/vguillemot/SHIP/reference/targetGstar.md)
  : Signed-correlation knowledge-based target.

## Data and helpers

- [`SHIP`](https://github.com/vguillemot/SHIP/reference/SHIP-package.md)
  [`SHIP-package`](https://github.com/vguillemot/SHIP/reference/SHIP-package.md)
  : SHIP provides shrinkage covariance estimation with user-selected
  targets. The available targets include diagonal, constant-correlation,
  and knowledge-based structures informed by functional gene groups.
- [`expl`](https://github.com/vguillemot/SHIP/reference/expl.md) : Small
  example extracted from a microarray data set.
- [`target.help()`](https://github.com/vguillemot/SHIP/reference/target.help.md)
  : Transform a list of Pathway IDs into a binary matrix.
- [`check.path()`](https://github.com/vguillemot/SHIP/reference/check.path.md)
  : Check if two genes belong to any common KEGG pathway.
