#' Log-likelihood for mapping populations
#' 
#' Log-likelihood for mapping populations
#'
#' The argument \code{counts} can be constructed using the \code{table} function for two markers. Genotype coding must represent dosage of a founder haplotype. For BC populations, possible allele dosages are 0,1. For DH and RIL pops, it is 0,2. For F2 and S1 pops, it is 0,1,2. S1r is an S1 population with the 1 alleles in repulsion phase.
#'
#' @param r recombination frequency
#' @param counts 3x3 contingency table for haplotype dosages 0,1,2
#' @param pop.type One of the following: "DH","BC","F2","S1r","RIL.self","RIL.sib"
#' 
#' @return log-likelihood
#' @export

LL <- function(r,counts,pop.type) {
  
  stopifnot(pop.type %in% c("DH","BC","F2","S1r","RIL.self","RIL.sib"))

  rownames(counts) <- colnames(counts) <- 0:2
  if (pop.type %in% c("DH","RIL.self","RIL.sib")) {
    R <- counts["0","2"] + counts["2","0"]  #N02 + N20
    NR <- counts["0","0"] + counts["2","2"] #N00 + N22
  }
  if (pop.type=="BC") {
    R <- counts["0","1"] + counts["1","0"]  #N02 + N20
    NR <- counts["0","0"] + counts["1","1"] #N00 + N22
  }
  if (pop.type=="F2") {
    N11 <- counts["1","1"]  #N11
    R.R <- counts["0","2"] + counts["2","0"]  #N02 + N20
    R.NR <- counts["0","1"] + counts["1","0"] + counts["1","2"] + counts["2","1"]  #N01 + N10 + N12 + N21
    NR.NR <- counts["0","0"] + counts["2","2"] #N00 + N22
  }
  if (pop.type=="S1r") {
    N11 <- counts["1","1"]  #N11
    NR.NR <- counts["0","2"] + counts["2","0"]  #N02 + N20
    R.NR <- counts["0","1"] + counts["1","0"] + counts["1","2"] + counts["2","1"]  #N01 + N10 + N12 + N21
    R.R <- counts["0","0"] + counts["2","2"] #N00 + N22
    pop.type <- "F2" #LL is the same in terms of R and NR
  }
  
  f1 <- function(r) {
    if (r > 0 & r <= 0.5) {
      s <- 1-r
      return(switch(pop.type,
             BC=R*log10(r)+NR*log10(s),
             DH=R*log10(r)+NR*log10(s),
             F2=2*R.R*log10(r)+2*NR.NR*log10(s)+R.NR*log10(r*s)+N11*log10(r^2+s^2),
             RIL.self=R*log10(r)-(R+NR)*log10(1+2*r),
             RIL.sib=R*log10(r)+NR*log10(1+2*r)-(R+NR)*log10(1+6*r)))
    } else {
      return(NA)
    }
  }
  return(apply(array(r),1,f1))
}
