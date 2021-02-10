#' Max Likelihood Estimation of Linkage
#' 
#' Max Likelihood Estimation of Linkage
#'
#' Can be used to estimate either the LOD score or recombination frequency, depending on the value of \code{LOD}. Genotype coding must represent dosage of a founder haplotype. For BC populations, possible allele dosages are {0,1}. For DH pops, it is {0,2}. For F2 pops, it is {0,1,2}.
#'
#' @param geno Matrix of haplotype dosages (markers x indiv)
#' @param pop.type One of the following: "DH","BC","F2"
#' @param LOD Logical, whether to return LOD (TRUE) or recomb freq (FALSE)
#' @param n.core For parallel execution on multiple cores
#' 
#' @return Matrix with RF or LOD
#' @export
#' @importFrom parallel makeCluster stopCluster parRapply clusterExport
#' @importFrom stats optimize

MLEL <- function(geno,pop.type,LOD,n.core=1) {
  
  m <- nrow(geno)
  tmp <- expand.grid(col=1:m,row=1:m)
  tmp <- tmp[tmp$row > tmp$col,]  #only need lower triangular
  
  cl <- makeCluster(n.core)
  clusterExport(cl=cl,varlist="LL")
  
  f1 <- function(x,geno,pop.type,LOD=FALSE) {
    counts <- table(factor(geno[x[1],],levels=0:2),factor(geno[x[2],],levels=0:2))
    ans <- optimize(f=LL,interval=c(0,0.5),counts=counts,pop.type=pop.type,maximum=TRUE)
    if (LOD) {
      tmp <- ans$objective-LL(r=0.5,counts=counts,pop.type=pop.type)
      return(max(tmp,0))
    } else {
      return(ans$maximum)
    }
  }
  
  ans <- parRapply(cl,tmp,f1,geno=geno,pop.type=pop.type,LOD=LOD)
  stopCluster(cl)
  
  outmat <- matrix(0,nrow=m,ncol=m)
  outmat[cbind(tmp[,1],tmp[,2])] <- ans
  outmat <- outmat + t(outmat)  #symmetric
  if (LOD) {
    diag(outmat) <- NA
  } else {
    diag(outmat) <- 1
  }
  colnames(outmat) <- rownames(outmat) <- rownames(geno)
  return(outmat)
}
