#' Phase S1 parent and reconstruct progeny in terms of parental haplotypes
#' 
#' Phase S1 parent and reconstruct progeny in terms of parental haplotypes
#' 
#' It is assumed that only segregating markers are present. Progeny reconstruction occurs using an HMM with a uniform transition probability matrix, based on an average recombination frequency \code{r}, and a uniform model for the genotype \code{error}.
#' 
#' @param geno ordered genotype matrix (markers x indiv) for one chromosome
#' @param r average recombination frequency to use for the HMM
#' @param error average genotype error to use for the HMM
#' 
#' @return List containing
#' \describe{
#' \item{parent}{two column matrix (rows = markers) with the haplotypes for the parent}
#' \item{progeny}{matrix with progeny reconstructed based on dosage of the second parental haplotype}
#' }
#' @export
#' @importFrom HMM initHMM viterbi

S1_haplo <- function(geno,r,error) {
  
  ans <- MLEL(geno=geno,pop.type="S1",LOD=FALSE,n.core=1,adjacent=TRUE)
  m <- nrow(geno)
  hap1 <- integer(m)
  for (i in 2:m) {
    if (ans$phase[i]=="c") {
      hap1[i] <- hap1[i-1]
    } else {
      hap1[i] <- abs(1-hap1[i-1])
    }
  }

  hap2 <- abs(1-hap1)
  parent.geno <- cbind(Hap1=hap1,Hap2=hap2)
  rownames(parent.geno) <- rownames(geno)
  
  #Recode based on dosage of haplotype 2
  geno2 <- t(apply(cbind(hap2,geno),1,function(u){
                                      if (u[1]==0) {
                                        return(abs(2-u[-1]))  
                                      } else {
                                        return(u[-1])
                                      }}))
                
  #Apply HMM
  s <- 1-r
  T.mat <- rbind(c(s^2,2*r*s,r^2),
              c(r*s,s^2+r^2,r*s),
              c(r^2,2*r*s,s^2))
  E.mat <- rbind(c(1-error,error/2,error/2),
                 c(error/2,1-error,error/2),
                 c(error/2,error/2,1-error))
  hmm <- initHMM(States=as.character(0:2),Symbols=as.character(0:2),
                 startProbs=c(1,2,1)/4,
                 transProbs=T.mat,
                 emissionProbs=E.mat)
  
  n <- ncol(geno2) #number of progeny
  geno.hmm <- geno2 #initialize to get same dimensions and names
  for (i in 1:n) {
    input <- as.character(geno2[,i])
    input[is.na(input)] <- "1" #this software does not handle missing data
    out <- viterbi(hmm,observation=input)  #returns ML solution
    geno.hmm[,i] <- as.integer(out) #convert from character back to integer
  }
  
  return(list(parent=parent.geno,progeny=geno.hmm))
}
