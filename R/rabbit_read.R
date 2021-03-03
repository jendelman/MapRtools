#' Parse output from RABBIT MagicReconstruct
#' 
#' Parse output from RABBIT MagicReconstruct
#' 
#' Two different file formats can be created. The \code{ML.file} contains the most likely (i.e., posterior maximum) genotype for each individual at each marker. The \code{diaQTL.file} contains the full distribution of genotype probabilities in the format required by the diaQTL R package (\code{diaQTL.file}). The default value for each filename is NULL, which generates no file.
#' 
#' @param rabbit.file name of RABBIT output file
#' @param ML.file name of most likely genotype file to create
#' @param diaQTL.file name of diaQTL genotype file to create
#' 
#' @return data frame defining the genotypes
#' @export
#' @importFrom utils write.csv

rabbit_read <- function(rabbit.file,ML.file=NULL,diaQTL.file=NULL) {

  stopifnot(!is.null(ML.file) | !is.null(diaQTL.file))
  con <- file(rabbit.file,"r",)
  suppressWarnings(temp <- readLines(con))
  close(con)
  
  ix <- grep("magicReconstruct",temp)

  k <- min(grep("genotype",temp[ix],ignore.case=F))
  x <- strsplit(temp[c((ix[k]+2):(ix[k+1]-1))],split=",")
  ng <- length(x)
  geno.code <- data.frame(code=1:ng,genotype=sapply(x,"[",3))

  k <- grep("genoprob",temp[ix],ignore.case=F)
  x <- strsplit(temp[c((ix[k]+1):(ix[k+1]-1))],split=",")
  map <- data.frame(marker=x[[1]][-1],chrom=x[[2]][-1],position=round(as.numeric(x[[3]][-1]),2),stringsAsFactors = F)

  y <- sapply(x[-(1:3)],function(x){round(as.numeric(x[-1]),3)})
  m <- nrow(map)
  n <- ncol(y)/ng

  w <- strsplit(sapply(x[-(1:3)],"[",1),split="_genotype",fixed=T)
  out1 <- matrix(0,nrow=m,ncol=n)
  colnames(out1) <- unique(sapply(w,"[",1))
  for (i in 1:m) {
    z <- split(y[i,],rep(1:n,each=ng))
    out1[i,] <- sapply(z,which.max)
  #out[i,] <- sapply(z,function(u){ix <- which(u > 0)
   #                     paste(paste((1:ng)[ix],collapse="|"),paste(u[ix],collapse="|"),sep="=>")})
  }

  if (!is.null(ML.file)) {
    write.csv(cbind(map,out1),file=ML.file,row.names=F)
  }
#write.csv(data.frame(id=colnames(out),parent1="F1",parent2="F1"),file="Week5/diaQTL_F2ped.csv",row.names=F)
  return(geno.code)
}