#' Creating a covariance target, optionally by using information from KEGG
#' pathways.
#' 
#' The function build.target() is a wrapper function to build the various types
#' of covariance targets (D,F,G,Gpos,Gstar,cor).
#' 
#' 
#' @param x An \eqn{n \times p}{n x p} matrix.
#' @param genegroups List of the groups each gene belongs to: each entry of the
#' list is dedicated to a gene (identified the same way as in \code{x}). Each
#' item of the list is thus a vector of pathway IDs.
#' @param type Character string specifying the wished target: "D" for a
#' diagonal target, "cor" for a correlation target, "G", "Gpos" and "Gstar" for
#' a G-type target (see Jelizarow et al, 2010) and "F" for a F-target.
#' @return A \eqn{p \times p}{p x p} target covariance matrix of a certain
#' type.
#' @author Vincent Guillemot
#' @seealso \code{\link{targetCor}}, \code{\link{targetD}},
#' \code{\link{targetF}}, \code{\link{targetG}}, \code{\link{targetGpos}},
#' \code{\link{targetGstar}},.
#' @references M. Jelizarow, V. Guillemot, A. Tenenhaus, K. Strimmer, A.-L.
#' Boulesteix, 2010.  Over-optimism in bioinformatics: an illustration.
#' Bioinformatics. Accepted.
#' @keywords methods
#' @examples
#' 
#' # Simulate dataset
#' x <- matrix(rnorm(20*30),20,30)
#' # Try different targets
#' build.target(x,type="D")
#' 
#' @export
build.target <- function(x, genegroups = NULL, type) {

  targetFun <- switch(type, cor = targetCor,
                            D = targetD,
                            F = targetF,
                            G = targetG,
                            Gpos = targetGpos,
                            Gstar = targetGstar)
  res <- targetFun(x,genegroups)
}
#' Computation of the target Cor.
#' 
#' The \eqn{p \times p}{p x p} target Cor is computed from the \eqn{n \times
#' p}{n x p} data matrix. It it a modified version of target G. In particular,
#' it tests the correlations (with a significance level of 0.05) and sets the
#' non-significant correlations to zero before the mean correlation
#' \eqn{\bar{r}}{r} is computed.
#' 
#' 
#' @param x A \eqn{n \times p} data matrix.
#' @param genegroups A list of genes obtained using the database KEGG, where
#' each entry itself is a list of pathway names this genes belongs to. If a
#' gene does not belong to any gene functional group, the entry is NA.
#' @return A \eqn{p \times p}{p x p} matrix.
#' @author Monika Jelizarow and Vincent Guillemot
#' @seealso \code{\link{targetCor}}, \code{\link{targetF}},
#' \code{\link{targetG}}, \code{\link{targetGstar}}, \code{\link{targetGpos}}.
#' @references J. Schaefer and K. Strimmer, 2005. A shrinkage approach to
#' large-scale covariance matrix estimation and implications for functional
#' genomics.  Statist. Appl. Genet. Mol. Biol. 4:32.
#' @keywords methods multivariate
#' @examples
#' 
#' # A short example on a toy dataset
#' # require(SHIP)
#' data(expl)
#' attach(expl)
#' tar <- targetCor(x,genegroups)
#' which(tar[upper.tri(tar)]!=0) # not many non zero coefficients !
#' @importFrom stats sd cov cor cor.test
#' @export
targetCor <- function(x,genegroups) {
    
    T1 <- target.help(genegroups)
    T2 <- matrix(nrow=length(genegroups),ncol=length(genegroups),data=0)     
    
    for (i in 2:length(genegroups)) {
      for (j in 1:(i-1)) if (T1[i,j]==1) T1[i,j] <- ifelse(stats::cor.test(x[,i],x[,j])$p.value < 0.05,1,0)
    }
    
    corm <- T1 * cor(x)     
    diag(corm) <- 0  
    
    cora <- ifelse(sum(corm!=0)==0,0,sum(colSums(corm))/sum(corm!=0))
    
    for (i in 1:length(genegroups)) {
      for (j in 1:i) {
        if (i!=j & T1[i,j]==1) T2[i,j] <- cora*(sd(x[,i])*sd(x[,j]))
        if (i==j)T2[i,j] <- cov(x[,i],x[,j])
        T2[j,i]<-T2[i,j]
      }
    }
    T2
  }
#' Computation of the diagonal target D ('diagonal, unequal variances').
#' 
#' The \eqn{p \times p}{n x p} diagonal target D is computed from the \eqn{n
#' \times p}{n x p} data matrix. It is defined as follows (\eqn{i,j =
#' 1,...,p}{i,j=1,...,p}): \deqn{t_{ij}=\left\{ }{tij = sii if i=j, 0
#' otherwise}\deqn{\begin {array} {ll} }{tij = sii if i=j, 0
#' otherwise}\deqn{s_{ii}\;&\mbox{if}\;i=j\\ }{tij = sii if i=j, 0
#' otherwise}\deqn{0\;&\mbox{if}\;i\neq j\\ }{tij = sii if i=j, 0
#' otherwise}\deqn{\end {array} }{tij = sii if i=j, 0
#' otherwise}\deqn{\right.}{tij = sii if i=j, 0 otherwise} where
#' \eqn{s_{ij}}{sij} denotes the entry of the unbiased covariance matrix in row
#' \eqn{i}{i}, column \eqn{j}{j}.
#' 
#' 
#' @param x A \eqn{n \times p}{n x p} data matrix.
#' @param genegroups The genegroups are not used for this target.
#' @return A \eqn{p \times p}{p x p} diagonal matrix.
#' @author Monika Jelizarow and Vincent Guillemot
#' @seealso \code{\link{targetCor}}, \code{\link{targetF}},
#' \code{\link{targetG}}, \code{\link{targetGstar}}, \code{\link{targetGpos}}.
#' @references J. Schaefer and K. Strimmer, 2005. A shrinkage approach to
#' large-scale covariance matrix estimation and implications for functional
#' genomics.  Statist. Appl. Genet. Mol. Biol. 4:32.
#' @keywords methods multivariate
#' @examples
#' 
#' x <- matrix(rnorm(10*30),10,30)
#' tar <- targetD(x,NULL)
#' 
#' @importFrom stats var
#' @export
targetD <- function(x,genegroups)  diag(apply(x,2,var))

#' Computation of target F ('constant correlation model').
#' 
#' The \eqn{p \times p}{p x p} target F is computed from the \eqn{n \times p}{n
#' x p} data matrix. It is defined as follows (\eqn{i,j = 1,...,p}{i,j =
#' 1,...,p}): \deqn{t_{ij}=\left\{ }{tij = sii if i=j ; r\sqrt{s_{ii}s_{jj}}
#' otherwise}\deqn{\begin {array} {ll} }{tij = sii if i=j ;
#' r\sqrt{s_{ii}s_{jj}} otherwise}\deqn{s_{ii}\;&\mbox{if}\;i=j\\ }{tij = sii
#' if i=j ; r\sqrt{s_{ii}s_{jj}}
#' otherwise}\deqn{\bar{r}\sqrt{s_{ii}s_{jj}}\;&\mbox{if}\;i\neq j\\ }{tij =
#' sii if i=j ; r\sqrt{s_{ii}s_{jj}} otherwise}\deqn{\end {array} }{tij = sii
#' if i=j ; r\sqrt{s_{ii}s_{jj}} otherwise}\deqn{\right.}{tij = sii if i=j ;
#' r\sqrt{s_{ii}s_{jj}} otherwise} where \eqn{\bar{r}}{r} is the average of
#' sample correlations and \eqn{s_{ij}}{sij} denotes the entry of the unbiased
#' covariance matrix in row \eqn{i}{i}, column \eqn{j}{j}.
#' 
#' 
#' @param x A \eqn{n \times p}{n x p} data matrix.
#' @param genegroups The genegroups are not used for this target.
#' @return A \eqn{p \times p}{p x p} matrix.
#' @author Monika Jelizarow and Vincent Guillemot
#' @seealso \code{\link{targetCor}}, \code{\link{targetF}},
#' \code{\link{targetG}}, \code{\link{targetGstar}}, \code{\link{targetGpos}}.
#' @references J. Schaefer and K. Strimmer, 2005. A shrinkage approach to
#' large-scale covariance matrix estimation and implications for functional
#' genomics.  Statist. Appl. Genet. Mol. Biol. 4:32.
#' @keywords methods multivariate
#' @examples
#' 
#' # A short example on a toy dataset
#' # require(SHIP)
#' data(expl)
#' attach(expl)
#' tar <- targetF(x,NULL)
#' which(tar[upper.tri(tar)]!=0) # many non zero coefficients !
#' 
#' @importFrom stats var
#' @export
targetF <-
  function(x,genegroups) {
    covm <- cov(x)
    corm <- cor(x)
    var.only <- covm/corm
    diag(corm) <- 0
    cora  <- sum(corm)/sum(corm!=0)
    T <- cora*var.only
    diag(T)<- diag(covm)
    T
  }

#' Computation of target G ('knowledge-based constant correlation model').
#' 
#' The \eqn{p \times p}{p x p} target G is computed from the \eqn{n \times p}{n
#' x p} data matrix. It is defined as follows (\eqn{i,j = 1,...,p}{i,j =
#' 1,...,p}): \deqn{t_{ij}=\left\{ }{sii if i=j ; r\sqrt{s_{ii}s_{jj}} if i~j ;
#' 0 otherwise }\deqn{\begin {array} {ll} }{sii if i=j ; r\sqrt{s_{ii}s_{jj}}
#' if i~j ; 0 otherwise }\deqn{s_{ii}\;&\mbox{if}\;i=j\\ }{sii if i=j ;
#' r\sqrt{s_{ii}s_{jj}} if i~j ; 0 otherwise
#' }\deqn{\bar{r}\sqrt{s_{ii}s_{jj}}\;&\mbox{if}\;i\neq j, i\sim j\\ }{sii if
#' i=j ; r\sqrt{s_{ii}s_{jj}} if i~j ; 0 otherwise }\deqn{0\;&\mbox{otherwise}
#' }{sii if i=j ; r\sqrt{s_{ii}s_{jj}} if i~j ; 0 otherwise }\deqn{\end{array}
#' }{sii if i=j ; r\sqrt{s_{ii}s_{jj}} if i~j ; 0 otherwise }\deqn{\right.}{sii
#' if i=j ; r\sqrt{s_{ii}s_{jj}} if i~j ; 0 otherwise } where \eqn{\bar{r}}{r}
#' is the average of sample correlations and \eqn{s_{ij}}{sij} denotes the
#' entry of the unbiased covariance matrix in row \eqn{i}{i}, column
#' \eqn{j}{j}. The notation \eqn{i\sim j}{i ~ j} means that genes \eqn{i}{i}
#' and \eqn{j}{j} are connected, i.e. genes \eqn{i}{i} and \eqn{j}{j} are in
#' the same gene functional group.
#' 
#' 
#' @param x A \eqn{n \times p}{n x p} data matrix.
#' @param genegroups A list of genes obtained using the database KEGG, where
#' each entry itself is a list of pathway names this genes belongs to. If a
#' gene does not belong to any gene functional group, the entry is NA.
#' @return A \eqn{p \times p}{p x p} matrix.
#' @author Monika Jelizarow and Vincent Guillemot
#' @seealso \code{\link{targetCor}}, \code{\link{targetF}},
#' \code{\link{targetG}}, \code{\link{targetGstar}}, \code{\link{targetGpos}}.
#' @references \itemize{ \item J. Schaefer and K. Strimmer, 2005. A shrinkage
#' approach to large-scale covariance matrix estimation and implications for
#' functional genomics.  Statist. Appl. Genet. Mol. Biol. 4:32.  \item M.
#' Jelizarow, V. Guillemot, A. Tenenhaus, K. Strimmer, A.-L. Boulesteix, 2010.
#' Over-optimism in bioinformatics: an illustration. Bioinformatics. Accepted.
#' }
#' @keywords methods multivariate
#' @examples
#' 
#' # A short example on a toy dataset
#' # require(SHIP)
#' data(expl)
#' attach(expl)
#' tar <- targetG(x,genegroups)
#' which(tar[upper.tri(tar)]!=0) # not many non zero coefficients !
#' 
#' @importFrom stats cov cor
#' @export
targetG <- function(x,genegroups) {
    T1   <- target.help(genegroups)
    T2   <- matrix(nrow=length(genegroups),ncol=length(genegroups),data=0)
    corm <- T1*cor(x)
    diag(corm)<-0
    cora <-sum(colSums(corm))/sum(corm!=0)
    for (i in 1:length(genegroups)) {
      for (j in 1:i) {
        if (i!=j & T1[i,j]==1) T2[i,j] <- cora*(sd(x[,i])*sd(x[,j]))
        if (i==j)T2[i,j] <- cov(x[,i],x[,j])
        T2[j,i]<-T2[i,j]
      }
    }
    T2
  }


#' Computation of the target Gpos.
#' 
#' The \eqn{p \times p}{p x p} target Gpos is computed from the \eqn{n \times
#' p}{n x p} data matrix. It it a modified version of target G. In particular,
#' it completely ignores negative correlations and computes the mean
#' correlation \eqn{\bar{r}}{r} using the positive ones only.
#' 
#' 
#' @param x A \eqn{n \times p}{n x p} data matrix.
#' @param genegroups A list of genes obtained using the database KEGG, where
#' each entry itself is a list of pathway names this genes belongs to. If a
#' gene does not belong to any gene functional group, the entry is NA.
#' @return A \eqn{p \times p}{p x p} matrix.
#' @author Monika Jelizarow and Vincent Guillemot
#' @seealso \code{\link{targetCor}}, \code{\link{targetF}},
#' \code{\link{targetG}}, \code{\link{targetGstar}}, \code{\link{targetGpos}}.
#' @references \itemize{ \item J. Schaefer and K. Strimmer, 2005. A shrinkage
#' approach to large-scale covariance matrix estimation and implications for
#' functional genomics.  Statist. Appl. Genet. Mol. Biol. 4:32. \item M.
#' Jelizarow, V. Guillemot, A. Tenenhaus, K. Strimmer, A.-L. Boulesteix, 2010.
#' Over-optimism in bioinformatics: an illustration. Bioinformatics. Accepted.
#' }
#' @keywords methods multivariate
#' @examples
#' 
#' # A short example on a toy dataset
#' # require(SHIP)
#' data(expl)
#' attach(expl)
#' tar <- targetGpos(x,genegroups)
#' which(tar[upper.tri(tar)]!=0) # not many non zero coefficients !
#' 
#' @importFrom stats sd cov cor
#' @export
targetGpos <- function(x,genegroups) {
    T1   <- target.help(genegroups)
    T2   <- matrix(nrow=length(genegroups),ncol=length(genegroups),data=0)
    corm <- T1*cor(x)
    diag(corm)<-0
    sum.pos <- sum(corm[T1==1 & corm > 0])
    cora.pos <- ifelse(sum(corm>0)==0,0,sum.pos/sum(corm>0))
    
    for (i in 1:length(genegroups)) {
      for (j in 1:i) {
        if (i!=j & T1[i,j]==1) T2[i,j] <- cora.pos*(sd(x[,i])*sd(x[,j]))
        if (i==j)T2[i,j] <- cov(x[,i],x[,j])
        T2[j,i]<-T2[i,j]
      }
    }
    T2
    
  }


#' Computation of the target Gstar.
#' 
#' The \eqn{p \times p}{p x p} target Gstar is computed from the \eqn{n \times
#' p}{n x p} data matrix. It it a modified version of target G. In particular,
#' it involves two parameters for the correlation (a positive and a negative
#' one) instead of the single parameter \eqn{\bar{r}}{r} in order to account
#' for negatively correlated genes within the same pathway
#' 
#' 
#' @param x A \eqn{n \times p} data matrix.
#' @param genegroups A list of genes obtained using the database KEGG, where
#' each entry itself is a list of pathway names this genes belongs to. If a
#' gene does not belong to any gene functional group, the entry is NA.
#' @return A \eqn{p \times p}{p x p} matrix.
#' @author Monika Jelizarow and Vincent Guillemot
#' @seealso \code{\link{targetCor}}, \code{\link{targetF}},
#' \code{\link{targetG}}, \code{\link{targetGstar}}, \code{\link{targetGpos}}.
#' @references \itemize{ \item J. Schaefer and K. Strimmer, 2005. A shrinkage
#' approach to large-scale covariance matrix estimation and implications for
#' functional genomics.  Statist. Appl. Genet. Mol. Biol. 4:32.  \item M.
#' Jelizarow, V. Guillemot, A. Tenenhaus, K. Strimmer, A.-L. Boulesteix, 2010.
#' Over-optimism in bioinformatics: an illustration. Bioinformatics. Accepted.
#' }
#' @keywords methods multivariate
#' @examples
#' 
#' # A short example on a toy dataset
#' # require(SHIP)
#' data(expl)
#' attach(expl)
#' tar <- targetGstar(x,genegroups)
#' which(tar[upper.tri(tar)]!=0) # not many non zero coefficients !
#' 
#' @importFrom stats sd cov cor
#' @export
targetGstar <- function(x,genegroups) {
    T1   <- target.help(genegroups)
    T2   <- matrix(nrow=length(genegroups),ncol=length(genegroups),data=0)
    corm <- T1*cor(x)
    diag(corm)<-0
    sum.pos <- sum(corm[T1==1 & corm > 0])
    sum.neg <- sum(corm[T1==1 & corm < 0])
    
    cora.pos <- sum.pos/sum(corm>0)
    cora.neg <- sum.neg/sum(corm<0)
    
    for (i in 1:length(genegroups)) {
      for (j in 1:i) {
        if (i!=j & T1[i,j]==1 & corm[i,j] > 0) T2[i,j] <- cora.pos*(sd(x[,i])*sd(x[,j]))
        if (i!=j & T1[i,j]==1 & corm[i,j] < 0) T2[i,j] <- cora.neg*(sd(x[,i])*sd(x[,j]))
        if (i==j)T2[i,j] <- cov(x[,i],x[,j])
        T2[j,i]<-T2[i,j]
      }
    }
    T2
  }

