#' Positive-correlation knowledge-based target.
#' 
#' This is a version of [targetG()] that estimates the common correlation from
#' positive linked correlations only. Negative linked correlations contribute
#' zero to the off-diagonal target.
#' \deqn{T_{ij} = \bar{r}_+\sqrt{s_{ii}s_{jj}}}{T_ij = r_plus * sqrt(s_ii * s_jj)}
#' for linked pairs with positive correlation, where \eqn{\bar{r}_+}{r_plus}
#' is their mean positive correlation.
#' 
#' 
#' @param x A \eqn{n \times p}{n x p} data matrix.
#' @param genegroups A list of genes obtained using the database KEGG, where
#' each entry itself is a list of pathway names this genes belongs to. If a
#' gene does not belong to any gene functional group, the entry is NA.
#' @return A \eqn{p \times p}{p x p} matrix.
#' @author Monika Jelizarow and Vincent Guillemot
#' @family covariance targets
#' @keywords methods multivariate
#' @examples
#' # require(SHIP)
#' data(expl)
#' target_matrix <- targetGpos(expl$x, expl$genegroups)
#' which(target_matrix[upper.tri(target_matrix)] != 0) # not many non zero coefficients !
#' 
#' @importFrom stats sd cov cor
#' @export
targetGpos <- function(x, genegroups) {
  validate_data(x)
  validate_gene_groups(genegroups, ncol(x))
  group_matrix <- target.help(genegroups)
  covariance_matrix <- cov(x)
  correlation_matrix <- group_matrix * cor(x)
  diag(correlation_matrix) <- 0
  positive_links <- group_matrix == 1 & correlation_matrix > 0
  positive_correlations <- correlation_matrix[positive_links & upper.tri(correlation_matrix)]
  positive_mean_correlation <- if (length(positive_correlations) == 0L) 0 else mean(positive_correlations)
  target_matrix <- positive_mean_correlation *
    tcrossprod(sqrt(diag(covariance_matrix))) * positive_links
  diag(target_matrix) <- diag(covariance_matrix)
  target_matrix
}
