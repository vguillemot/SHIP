#' Transform a list of Pathway IDs into a binary matrix.
#' 
#' This function is an internal function used to transform the given list of
#' \eqn{p}{p} (one vector of pathway IDs per gene) groups into a binary matrix.
#' 
#' 
#' @param genes List of \eqn{p}{p} items. Each item is the vector of Pathway
#' IDs a gene belongs to.
#' @return A \eqn{p \times p}{p x p} binary matrix: the coefficient (i,j) is 1
#' if genes i and j belong to a common pathway and 0 otherwise. This is an
#' internal function called by for example \code{\link{targetG}}.
#' @author Monika Jelizarow and Vincent Guillemot
#' @seealso
#' \code{\link{targetF}},\code{\link{targetG}},\code{\link{targetGpos}},
#' \code{\link{targetGstar}}.
#' @keywords internal
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
  T <- diag(1,length(genes))                              
  for (i in 2:length(genes)) {                            
   for(j in 1:(i-1)) {                                    
    T[j,i] <- T[i,j] <- check.path(genes[[i]],genes[[j]]) 
   }
  }
  return(T)
}

