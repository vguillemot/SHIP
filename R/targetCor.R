#' Significance-filtered knowledge-based target.
#' 
#' The \eqn{p \times p}{p x p} target Cor is computed from the \eqn{n \times
#' p}{n x p} data matrix. It it a modified version of target G. In particular,
#' it tests the correlations (with a significance level of 0.05) and sets the
#' non-significant correlations to zero before the mean correlation
#' \eqn{\bar{r}}{r} is computed.
#'
#' For retained linked pairs, the off-diagonal target is
#' \deqn{T_{ij} = \bar{r}_{sig}\sqrt{s_{ii}s_{jj}}}{T_ij = r_sig * sqrt(s_ii * s_jj)}
#' where \eqn{\bar{r}_{sig}}{r_sig} is the mean correlation among links
#' whose correlation test p-value is below 0.05.
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
#' 
#' # A short example on a toy dataset
#' # require(SHIP)
#' data(expl)
#' target_matrix <- targetCor(expl$x, expl$genegroups)
#' which(target_matrix[upper.tri(target_matrix)] != 0) # not many non zero coefficients !
#' @importFrom stats sd cov cor cor.test
#' @export
targetCor <- function(x, genegroups) {
  validate_data(x)
  validate_gene_groups(genegroups, ncol(x))

  group_matrix <- target.help(genegroups)
  covariance_matrix <- cov(x)
  correlation_matrix <- cor(x)
  if (length(genegroups) > 1L) {
    for (i in 2:length(genegroups)) {
      for (j in 1:(i - 1L)) {
        if (group_matrix[i, j] == 1) {
          significant <- stats::cor.test(x[, i], x[, j])$p.value < 0.05
          group_matrix[i, j] <- as.integer(significant)
          group_matrix[j, i] <- group_matrix[i, j]
        }
      }
    }
  }

  filtered_correlations <- correlation_matrix[group_matrix == 1 & upper.tri(correlation_matrix)]
  mean_correlation <- if (length(filtered_correlations) == 0L) 0 else mean(filtered_correlations)
  target_matrix <- mean_correlation *
    tcrossprod(sqrt(diag(covariance_matrix))) * group_matrix
  diag(target_matrix) <- diag(covariance_matrix)
  target_matrix
}
