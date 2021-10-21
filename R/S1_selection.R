#' Signatures of selection in S1 populations
#' 
#' Signatures of selection in S1 populations
#' 
#' The null hypothesis is no selection, in which case the expected frequency of genotypes is (AA = 1/4, AB = 1/2, BB = 1/4). Two alternate hypotheses, or "modes", are tested for gametic selection: 1.selection in one sex, 2.selection in both sexes. Two modes of zygotic selection are also tested: 1.selection against one homozygote, 2.selection against both homozygotes. In addition to these 4 models with pure gametic or zygotic selection, 4 models with both forms of selection are considered (2 gametic modes x 2 zygotic modes = 4 models). The selection coefficients and -log10(p) value are returned for the model with the lowest AIC. The selection coefficients represent the sum of the absolute differences between the observed and expected frequencies. Positive values correspond to selection against A or AA, negative values for selection against B or BB. For zygotic mode 2, positive values represent selection against homozygotes, while negative values indicate selection against heterozygotes.  P-values are computed based on the likelihood ratio test; in other words, the change in deviance is assumed to be chi-squared distributed under the null hypothesis.
#' 
#' @param x Vector of progeny counts for the three possible genotypes (AA, AB, BB)
#' 
#' @return List containing
#' \describe{
#' \item{s}{selection coefficients}
#' \item{model}{the gametic and zygotic selection modes (NA,1,2) for the best model}
#' \item{score}{-log10(p) value}
#' }
#' 
#' @importFrom stats optimize pchisq
#' @importFrom utils capture.output
#' @import optimx
#' @export
#' 

S1_selection <- function(x) {
  
  freq <- function(s,mode) {

    gamete.mode <- mode[1]
    zygote.mode <- mode[2]
    
    k <- which(is.na(mode))
    if (length(k) > 0) {
      #one selection coefficient
      if (k==1) {
        sg <- 0
        sz <- s
        gamete.mode <- 1
      } else {
        sz <- 0
        sg <- s
        zygote.mode <- 1
      }
    } else {
      sg <- s[1]
      sz <- s[2]
    }
    
    if (gamete.mode==1) {
      p1 <- 1/2 - sg/2
      p2 <- 1/2
    } else {
      p1 <- 1/2 - sg/2
      p2 <- 1/2 - sg/2
    }
    q1 <- 1 - p1
    q2 <- 1 - p2
    
    f <- function(s) {
      ifelse(s > 0,s/2,s/4)
    }
    
    #Zygote frequencies
    P <- p1*p2
    H <- p1*q2 + q1*p2
    Q <- q1*q2
    
    if (zygote.mode==1) {
      P <- P - f(sz)
      H <- H + abs(sz)/4
      Q <- Q - f(-sz)
    } else {
      P <- P - sz/4
      H <- H + sz/2
      Q <- Q - sz/4
    }
    return(c(P,H,Q))
  }
  
  LL <- function(s,x,mode) {
    p <- freq(s,mode)
    return(sum(x*log(p)))
  }
  
  cases <- expand.grid(1:2,1:2)
  colnames(cases) <- c("gamete","zygote")
  cases <- split(rbind(cases,c(1,NA),c(2,NA),c(NA,1),c(NA,2)),1:8)
  suppressWarnings(ans <- lapply(cases,function(mode){
    mode <- as.integer(mode)
    if (any(is.na(mode))) {
      capture.output(ans <- optimx(par=0,fn=LL,gr=NULL,lower=-1,upper=1,
             x=x,mode=mode,control=list(maximize=TRUE),method=c("L-BFGS-B","nlminb","Rcgmin","Rvmmin")))
      k <- which.max(ans$value)
      list(ans[k,1],ans[k,2])
    } else {
      capture.output(ans <- optimx(par=c(0,0),fn=LL,gr=NULL,lower=c(-1,-1),upper=c(1,1),
             x=x,mode=mode,control=list(maximize=TRUE),method=c("L-BFGS-B","nlminb","Rcgmin","Rvmmin")))
      k <- which.max(ans$value)
      list(as.numeric(ans[k,1:2]),ans[k,3])
    }
  }))
  
  AIC <- 2*c(2,2,2,2,1,1,1,1) - 2*sapply(ans,"[[",2)
  
  LL0 <- LL(s=c(0,0),x=x,mode=c(1,1))
  k <- which.min(AIC)
  LL1 <- as.numeric(ans[[k]][2])
  soln <- ans[[k]][[1]]
  score <- -log10(pchisq(q=2*(LL1-LL0),df=length(soln),lower.tail=F))
  model <- as.numeric(cases[[k]])
  names(model) <- c("gametic","zygotic")
  return(list(s=soln,model=model,score=score))
}
