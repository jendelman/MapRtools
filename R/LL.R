#' Log-likelihood for inbred line-derived mapping populations
#' 
#' Log-likelihood for inbred line-derived mapping populations
#'
#' The argument \code{counts} can be constructed using the \code{table} function for two markers. Genotype coding must represent dosage of a founder haplotype. For BC populations, possible allele dosages are {0,1}. For DH pops, it is {0,2}. For F2 pops, it is {0,1,2}.
#'
#' @param r recombination frequency
#' @param counts 3x3 contingency table for haplotype dosages 0,1,2
#' @param pop.type One of the following: "DH","BC","F2"
#' 
#' @return log-likelihood
#' @export

LL <- function(r,counts,pop.type) {
  
  rownames(counts) <- colnames(counts) <- 0:2
  if (pop.type %in% "DH") {
    N.r <- counts["0","2"] + counts["2","0"]  #N02 + N20
    N.nr <- counts["0","0"] + counts["2","2"] #N00 + N22
  }
  if (pop.type=="BC") {
    N.r <- counts["0","1"] + counts["1","0"]  #N02 + N20
    N.nr <- counts["0","0"] + counts["1","1"] #N00 + N22
  }
  if (pop.type=="F2") {
    N11 <- counts["1","1"]  #N11
    N.r.r <- counts["0","2"] + counts["2","0"]  #N02 + N20
    N.r.nr <- counts["0","1"] + counts["1","0"] + counts["1","2"] + counts["2","1"]  #N01 + N10 + N12 + N21
    N.nr.nr <- counts["0","0"] + counts["2","2"] #N00 + N22
  }
  
  f1 <- function(r) {
    if (r > 0 & r <= 0.5) {
      s <- 1-r
      if (pop.type=="F2") {
        return(2*N.r.r*log10(r)+2*N.nr.nr*log10(s)+N.r.nr*log10(r*s)+N11*log10(r^2+s^2))
      } else {
        return(N.r*log10(r)+N.nr*log10(s))
      }
    } else {
      return(NA)
    }
  }
  return(apply(array(r),1,f1))
}
