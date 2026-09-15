#' Knowledge-based constant-correlation target.
#' 
#' This target applies the constant-correlation model only to pairs of genes
#' that share a functional group. Let \eqn{C = (c_{ij})}{C = (c_ij)} be the
#' binary group matrix and let \eqn{\bar{r}}{r_bar} be the mean correlation
#' among linked pairs. Then
#' \deqn{t_{ij} = 
#'   \begin{cases} 
#'     s_{ii} & \text{ if } i=j\\
#'     \bar{r}\sqrt{s_{ii}s_{jj}} & \text{ if } i\neq j, i\sim j
#'   \end{cases}} 
#' where \eqn{i \sim j}{i ~ j} means that genes \eqn{i}{i} and \eqn{j}{j}
#' share a functional group. Unlinked off-diagonal entries are zero.
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
#' target_matrix <- targetG(expl$x, expl$genegroups)
#' which(target_matrix[upper.tri(target_matrix)] != 0) # not many non zero coefficients !
#' 
#' @importFrom stats cov cor
#' @export
targetG <- function(x, genegroups) {
  validate_data(x)
  validate_gene_groups(genegroups, ncol(x))
  group_matrix <- target.help(genegroups)
  covariance_matrix <- cov(x)
  correlation_matrix <- group_matrix * cor(x)
  diag(correlation_matrix) <- 0
  correlations <- correlation_matrix[group_matrix == 1 & upper.tri(correlation_matrix)]
  mean_correlation <- if (length(correlations) == 0L) 0 else mean(correlations)
  target_matrix <- mean_correlation *
    tcrossprod(sqrt(diag(covariance_matrix))) * group_matrix
  diag(target_matrix) <- diag(covariance_matrix)
  target_matrix
}
