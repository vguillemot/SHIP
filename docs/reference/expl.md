# Small example extracted from a microarray data set.

The microarray data set is the study on the prostate cancer by Singh et
al. The collection of the microarray is hgu95av2, and the gene groups
are thus given by the information in the hgu95av2.db Bioconductor
library (see Carlson et al.).

## Usage

``` r
data("expl")
```

## Format

The dataset is a list containing:

- a \\102 \times 100\\ matrix \\x\\ of 100 genes randomly chosen from
  the data set of Singh et al.,

- a list `genegroups` containing 100 vectors of KEGG pathway IDs (which
  each gene belongs to).

## Source

- M. Carlson, S. Falcon, H. Pages, N. Li. hgu95av2.db: Affymetrix Human
  Genome U95 Set annotation data (chip hgu95av2). R package version
  2.2.12.

- D. Singh, P. G. Febbo, K. Ross, D. G. Jackson, J. Manola, C. Ladd, P.
  Tamayo, A. A. Renshaw, A. V. D'Amico, J. P. Richie, E. S. Lander, M.
  Loda, P. W. Kantoff, T. R. Golub, W. R. Sellers, 2002. Gene expression
  correlates of clinical prostate cancer behavior. Cancer Cell,
  Department of Adult Oncology, Brigham and Women's Hospital, Harvard
  Medical School, Boston, MA 02115, USA., 1, 203-209.
