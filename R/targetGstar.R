#' Signed-correlation knowledge-based target.
#' 
#' The \eqn{p \times p}{p x p} target Gstar is computed from the \eqn{n \times
#' p}{n x p} data matrix. It it a modified version of target G. In particular,
#' it involves two parameters for the correlation (a positive and a negative
#' one) instead of the single parameter \eqn{\bar{r}}{r} in order to account
#' for negatively correlated genes within the same pathway.
#'
#' The target uses separate mean correlations for positive and negative linked
#' pairs:
#' \deqn{T_{ij} = \begin{cases}
#'   \bar{r}_+\sqrt{s_{ii}s_{jj}} & \text{if } r_{ij} > 0 \\
#'   \bar{r}_-\sqrt{s_{ii}s_{jj}} & \text{if } r_{ij} < 0
#' \end{cases}}{T_ij = r_plus * sqrt(s_ii * s_jj) if r_ij > 0, and r_minus * sqrt(s_ii * s_jj) if r_ij < 0.}
#' Here \eqn{\bar{r}_+}{r_plus} and \eqn{\bar{r}_-}{r_minus} are the mean
#' positive and negative linked correlations.
#' 
#' 
#' @param x A \eqn{n \times p} data matrix.
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
#' target_matrix <- targetGstar(expl$x, expl$genegroups)
#' which(target_matrix[upper.tri(target_matrix)] != 0) # not many non zero coefficients !
#' 
#' @importFrom stats sd cov cor
#' @export
targetGstar <- function(x, genegroups) {
  validate_data(x)
  validate_gene_groups(genegroups, ncol(x))
  group_matrix <- target.help(genegroups)
  covariance_matrix <- cov(x)
  correlation_matrix <- group_matrix * cor(x)
  diag(correlation_matrix) <- 0
  
  positive_links <- group_matrix == 1 & correlation_matrix > 0
  negative_links <- group_matrix == 1 & correlation_matrix < 0
  positive_correlations <- correlation_matrix[positive_links & upper.tri(correlation_matrix)]
  negative_correlations <- correlation_matrix[negative_links & upper.tri(correlation_matrix)]
  positive_mean_correlation <- if (length(positive_correlations) == 0L) 0 else mean(positive_correlations)
  negative_mean_correlation <- if (length(negative_correlations) == 0L) 0 else mean(negative_correlations)
  target_matrix <- tcrossprod(sqrt(diag(covariance_matrix))) *
    (positive_mean_correlation * positive_links +
      negative_mean_correlation * negative_links)
  diag(target_matrix) <- diag(covariance_matrix)
  target_matrix
}

