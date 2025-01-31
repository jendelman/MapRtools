#' Expected Genotype Quality
#' 
#' Expected Genotype Quality for Binomial Model
#' 
#' Expected GQ conditional on the true allele dosage (0,1,2,...ploidy)
#' 
#' @param DP read depth
#' @param error allelic error
#' @param prior numeric vector of length ploidy + 1
#' 
#' @return numeric vector of length ploidy + 1
#' @export
#' @importFrom stats dbinom

EGQ <- function(DP,error,prior) {

  likelihood <- function(n,eps,ploidy){
    ans <- NULL
    for (z in 0:ploidy) {
      psi <- z/ploidy
      ans <- rbind(ans,
                   data.frame(z=z,
                              k=0:n,
                              prob=dbinom(x=0:n,size=n,prob=psi*(1-eps)+(1-psi)*eps)))
    }
    return(ans)
  }

  posterior <- function(n,eps,ploidy,prior) {
    tmp <- likelihood(n,eps,ploidy)
    tmp$posterior <- tmp$prob*rep(prior,each=n+1)
    denom <- rep(tapply(tmp$posterior,tmp$k,sum),times=ploidy+1)
    tmp$posterior <- tmp$posterior/denom
    return(tmp)
  }
  ploidy <- length(prior) - 1
  tmp <- posterior(n=DP,eps=error,ploidy=ploidy,prior=prior)
  GQ <- -10*log10(tapply(tmp$posterior,tmp$k,function(z){1-max(z)}))
  tmp <- likelihood(n=DP,eps=error,ploidy=ploidy)
  #sum(GQ*tapply(tmp$prob*rep(prior,each=DP+1),tmp$k,sum))
  tapply(rep(GQ,ploidy+1)*tmp$prob,tmp$z,sum)
}

