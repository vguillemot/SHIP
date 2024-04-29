\name{build.target}
\alias{build.target}
\title{
Creating a covariance target, optionally by using information from KEGG pathways.
}
\description{
The function build.target() is a wrapper function to build the various types of
covariance targets (D,F,G,Gpos,Gstar,cor).
}
\usage{ 
build.target(x, genegroups = NULL, type)
}

\arguments{
  \item{x}{ An \eqn{n \times p}{n x p} matrix.
}
  \item{genegroups}{
List of the groups each gene belongs to: each entry of the list is dedicated to a 
gene (identified the same way as in \code{x}). Each item of the list is thus a 
vector of pathway IDs.
}
  \item{type}{Character string specifying the wished target: "D" for a diagonal target, 
  "cor" for a correlation target, "G", "Gpos" and "Gstar" for a G-type target 
  (see Jelizarow et al, 2010) and "F" for a F-target.
      
}
}
\value{
A \eqn{p \times p}{p x p} target covariance matrix of a certain type.
}
\references{
 M. Jelizarow, V. Guillemot, A. Tenenhaus, K. Strimmer, A.-L. Boulesteix, 2010. 
 Over-optimism in bioinformatics: an illustration. Bioinformatics. Accepted.
}
\author{
Vincent Guillemot
}
\seealso{
\code{\link{targetCor}}, \code{\link{targetD}}, \code{\link{targetF}}, 
\code{\link{targetG}}, \code{\link{targetGpos}}, \code{\link{targetGstar}},.
}
\examples{
# Simulate dataset
x <- matrix(rnorm(20*30),20,30)
# Try different targets
build.target(x,type="D")

}
\keyword{methods}

