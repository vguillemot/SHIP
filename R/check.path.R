#' Check if two genes belong to any common KEGG pathway.
#' 
#' Takes as arguments two vectors of IDs and test whether they have a common
#' ID.
#' 
#' 
#' @param p1 Vector of pathways that gene 1 belongs to.
#' @param p2 Vector of pathways that gene 2 belongs to.
#' @return Return 0 if the two genes don't belong to a common pathway, return 1
#' otherwise.  This is an internal function used by the function
#' \code{\link{target.help}}.
#' @author Monika Jelizarow and Vincent Guillemot
#' @seealso \code{\link{target.help}}
#' @keywords internal
#' @examples
#' 
#' g1 <- c("path1","path2","path3","path4")
#' g2 <- c("path5","path6","path3","path11")
#' g3 <- c("path10","path5","path12","path13")
#' check.path(g1, g2) # 1
#' check.path(g1, g3) # 0
#' @export
check.path <- function(p1,p2) {
  bool <- any(p1 %in% p2, na.rm = TRUE) + 0
  return(bool)
}

