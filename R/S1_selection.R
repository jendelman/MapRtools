#' Signatures of selection in S1 populations
#' 
#' Signatures of selection in S1 populations
#' 
#' The null hypothesis is no selection, in which case the expected frequency of genotypes is (AA = 1/4, AB = 1/2, BB = 1/4). Four alternate hypotheses are tested, corresponding to 1. gametic selection in one sex, 2. gametic selection in both sexes, 3. zygotic selection against one homozygote, 4. zygotic selection against both homozygotes. The selection coefficient (s) is estimated by maximum likelihood on the interval [-1,1], where positive values correspond to selection against A or AA and negative values are for selection against B or BB. For model 4, positive values represent selection against homozygotes, while negative values indicate selection against heterozygotes.  P-values are computed based on the likelihood ratio test; in other words, the change in deviance is assumed to be chi-squared (df = 1) distributed under the null hypothesis.
#' 
#' @param x Vector of progeny counts for the three possible genotypes (AA, AB, BB)
#' 
#' @return Data frame with columns "s" for the ML selection coefficient and "score" for the -log10(p) value.
#' 
#' @importFrom stats optimize pchisq
#' @export
#' 

S1_selection <- function(x) {
  
  model1 <- function(s) {
    return(c(1-s,2,1+s)/4)
  }
  model2 <- function(s) {
    return(c((1-s)^2,2*(1-s^2),(1+s)^2)/4)
  }
  f <- function(s) {
    ifelse(s > 0,s,s/3)
  }
  model3 <- function(s) {
    return(c(1-f(s),2*(1+abs(s)/3),1-f(-s))/4)
  }
  model4 <- function(s) {
    return(c(1-s,2*(1+s),1-s)/4)
  }
  
  LL <- function(s,x,model) {
    p <- do.call(paste0("model",model),args=list(s=s))
    return(sum(x*log(p)))
  }
  ans <- lapply(list(1,2,3,4),FUN=function(model){optimize(f=LL,interval=c(-0.9999,0.9999),x=x,model=model,maximum = T)})
  LL0 <- sum(x*log(model1(s=0)))
  LL1 <- sapply(ans,"[[",'objective')
  data.frame(s=sapply(ans,"[[",'maximum'),
             score=-log10(1-pchisq(q=2*(LL1-LL0),df=1,lower.tail=T)))
}
