#' Convert progeny genotypes to reflect parental haplotypes
#' 
#' Convert progeny genotypes to reflect parental haplotypes
#'
#' Converts progeny genotypes so that the "0" allele reflects the haplotype inherited from parent 1 and the "1" allele reflects the haplotype from parent 2.
#'
#' @param geno A matrix of genotypes (markers x individuals) coded "0/0", "1/1" or "0/1"
#' @param par1 chr indicating the column name of parent 1
#' @param par2 chr indicating the column name of parent 2

#' @return A matrix with parental lines removed, and genotypes recoded to reflect parental haplotypes
#' @export

recodeGenotypes <- function(geno, par1, par2){
  par1.idx <- which(colnames(geno) == par1)
  par2.idx <- which(colnames(geno) == par2)
  
  geno.converted <- matrix(nrow = nrow(geno), ncol = ncol(geno))
  colnames(geno.converted) <- colnames(geno)
  rownames(geno.converted) <- rownames(geno)
  
  mrk.idx <- 1
  
  for(m in 1:nrow(geno)){
    par1.geno <- geno[m,par1.idx]
    par2.geno <- geno[m,par2.idx]
    
    if(!is.na(par1.geno) & !is.na(par2.geno)){
      
      if(par1.geno == "0/0" && par2.geno == "1/1"){
        geno.converted[mrk.idx, ] <- geno[m,]
        rownames(geno.converted)[mrk.idx] <- rownames(geno)[m]
        mrk.idx <- mrk.idx + 1
      }
      else if(par1.geno == "1/1" && par2.geno == "0/0"){
        geno.converted[mrk.idx, ] <- chartr("01", "10", geno[m, ])
        rownames(geno.converted)[mrk.idx] <- rownames(geno)[m]
        mrk.idx <- mrk.idx + 1
      }
    }
  }
  geno.converted <- geno.converted[,-c(par1.idx, par2.idx)]
  geno.converted <- geno.converted[-c(mrk.idx:m),]
  return(geno.converted)
}