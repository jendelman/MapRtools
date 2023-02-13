#' Expected Genotype Quality
#' 
#' Expected Genotype Quality for Binomial Model
#' 
#' As defined in Matias et al. (2019), EGQ is the PHRED-scaled expected error of the genotype call, conditional on the true genotype. This function returns EGQ for the genotype most frequently miscalled, which is the balanced heterozygote (i.e., ploidy/2).
#' 
#' @param depth read count
#' @param error allelic error
#' @param ploidy ploidy
#' @param prior numeric vector of length ploidy + 1
#' 
#' @return numeric scalar
#' @export
#' @references Matias et al. (2019) Plant Genome 12:190002. https://doi.org/10.3835/plantgenome2019.01.0002

EGQ <- function(depth,error,ploidy,prior) {
  
  posterior <- function(k,N,error,ploidy,prior){
    p <- seq(0,1,by=1/ploidy)
    q <- 1-p
    z <- dbinom(x=k,size=N,prob=p*(1-error)+q*error)*prior
    names(z) <- 0:ploidy
    return(z/sum(z))
  }
  
  f1 <- function(N,error,ploidy,prior) {
    k <- seq(0,N,by=1)
    p <- apply(array(k),1,posterior,N=N,error=error,ploidy=ploidy,prior=prior)
    geno.mode <- apply(p,2,which.max)-1
    prob.error <- sum(dbinom(k,N,prob=0.5)[which(geno.mode!=ploidy/2)])
    return(-10*log10(prob.error))
  }
  min(apply(array(depth:(depth+10)),1,f1,error=error,ploidy=ploidy,prior=prior))
}
