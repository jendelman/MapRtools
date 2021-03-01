#' Max Likelihood Estimation of Linkage
#' 
#' Max Likelihood Estimation of Linkage
#'
#' Can be used to estimate either the LOD score or recombination frequency, depending on the value of \code{LOD}. Genotype coding must represent dosage of a founder haplotype. For BC populations, possible allele dosages are 0,1. For DH and RIL pops, it is 0,2. For F2 and S1 pops, it is 0,1,2. 
#'
#' @param geno Matrix of haplotype dosages (markers x indiv)
#' @param pop.type One of the following: "DH","BC","F2","S1","RIL.self","RIL.sib"
#' @param LOD Logical, whether to return LOD (TRUE) or recomb freq (FALSE)
#' @param n.core For parallel execution on multiple cores
#' @param adjacent Logical, should calculation be done for all pairs (FALSE) or adjacent (TRUE) markers
#' 
#' @return If \code{adjacent} is FALSE, a matrix of recombination frequencies or LOD scores; otherwise, a three-column data frame with marker, the LOD or r value, and the phase ("c","r") with the previous marker 
#' @export
#' @importFrom parallel makeCluster stopCluster parRapply clusterExport
#' @importFrom stats optimize

MLEL <- function(geno,pop.type,LOD,n.core=1,adjacent=FALSE) {
  
  f1 <- function(x,geno,pop.type,LOD=FALSE) {
    counts <- table(factor(geno[x[1],],levels=0:2),factor(geno[x[2],],levels=0:2))
    
    if (pop.type %in% c("F2","S1")) {
      ans <- optimize(f=LL,interval=c(0,0.5),counts=counts,pop.type="F2",maximum=TRUE)  
      phase <- 0 #coupling
      if (pop.type=="S1") {
        ans.r <- optimize(f=LL,interval=c(0,0.5),counts=counts,pop.type="S1r",maximum=TRUE) #repulsion
        #choose phase based on ML
        if (ans$objective < ans.r$objective) {
          ans <- ans.r
          phase <- 1 #repulsion
        } 
      } 
      if (LOD) {
        tmp <- ans$objective-LL(r=0.5,counts=counts,pop.type="F2")
        result <- max(tmp,0)
      } else {
        result <- ans$maximum
      }
    } else {
      phase <- 0  #coupling assumed
      if (pop.type %in% c("DH","RIL.self","RIL.sib")) {
        R <- counts["0","2"] + counts["2","0"]  #N02 + N20
        NR <- counts["0","0"] + counts["2","2"] #N00 + N22
      }
      if (pop.type=="BC") {
        R <- counts["0","1"] + counts["1","0"]  #N02 + N20
        NR <- counts["0","0"] + counts["1","1"] #N00 + N22
      }
      
      if (pop.type %in% c("BC","DH")) {
        rML <- R/(R+NR)
      }
      if (pop.type=="RIL.self") {
        rML <- R/2/NR
      }
      if (pop.type=="RIL.sib") {
        rML <- R/(4*NR-2*R)
      }
      rML <- min(rML,0.5)
      if (LOD) {
        if (rML==0) {
          tmp <- Inf 
        } else {
          tmp <- LL(r=rML,counts=counts,pop.type=pop.type) - LL(r=0.5,counts=counts,pop.type=pop.type)
        }
        result <- max(tmp,0)
      } else {
        result <- rML
      }
    }
    return(c(result,phase))
  }
  
  m <- nrow(geno)
  if (!adjacent) {
    tmp <- expand.grid(col=1:m,row=1:m)
    tmp <- tmp[tmp$row > tmp$col,]  #only need lower triangular
  } else {
    tmp <- cbind(1:(m-1),2:m)
  }
  
  cl <- makeCluster(n.core)
  clusterExport(cl=cl,varlist="LL")
  ans <- parRapply(cl,tmp,f1,geno=geno,pop.type=pop.type,LOD=LOD)
  stopCluster(cl)
  ans <- matrix(ans,byrow=TRUE,ncol=2)  #first column is result, second is phase
  
  if (!adjacent) {
    out <- matrix(0,nrow=m,ncol=m)
    out[cbind(tmp[,1],tmp[,2])] <- ans[,1]
    out <- out + t(out)  #symmetric
    if (LOD) {
      diag(out) <- Inf
    } else {
      diag(out) <- 0
    }
    colnames(out) <- rownames(out) <- rownames(geno)
    return(out)
  } else {
    out <- data.frame(marker=rownames(geno),
                      value=c(NA,ans[,1]),
                      phase=c(as.character(NA),ifelse(ans[,2]==0,"c","r")))
    return(out)
  }
}
