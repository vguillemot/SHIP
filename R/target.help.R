#' Transform a list of Pathway IDs into a binary matrix.
#' 
#' This function transforms a list of \eqn{p}{p} (one vector of pathway IDs per
#' gene) groups into a binary matrix.
#' 
#' 
#' @param genes List of \eqn{p}{p} items. Each item is the vector of Pathway
#' IDs a gene belongs to.
#' @return A \eqn{p \times p}{p x p} binary matrix: the coefficient (i,j) is 1
#' if genes i and j belong to a common pathway and 0 otherwise.
#' @author Monika Jelizarow and Vincent Guillemot
#' @seealso
#' \code{\link{targetF}},\code{\link{targetG}},\code{\link{targetGpos}},
#' \code{\link{targetGstar}}.
#' @keywords methods
#' @examples
#' 
#' g1 <- c("path1", "path2", "path3", "path4")
#' g2 <- c("path5", "path6", "path3", "path11")
#' g3 <- c("path10", "path5", "path12", "path13")
#' target.help(list(g1, g2, g3)) 
#' 
#' 
#' @export
target.help <- function(genes) {
  validate_gene_groups(genes, length(genes))
  group_matrix <- diag(1, length(genes))
  if (length(genes) < 2L) {
    return(group_matrix)
  }
  for (i in 2:length(genes)) {
    for (j in 1:(i - 1L)) {
      group_matrix[j, i] <- group_matrix[i, j] <- check.path(genes[[i]], genes[[j]])
    }
  }
  group_matrix
}

