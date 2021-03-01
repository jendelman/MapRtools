#' Parse output from RABBIT MagicReconstruct
#' 
#' Parse output from RABBIT MagicReconstruct
#' 
#' Generates output files with the maximum likelihood (ML) genotype (\code{out1file}) and genotype probabilities for the diaQTL package (\code{out2file}).
#' 
#' @param filename name of RABBIT output file
#' @param outfile1 name of ML genotype file to create
#' @param outfile2 name of diaQTL genotype file to create 
#' 
#' @export

rabbit_read <- function(filename,outfile1="MLgeno.csv",outfile2=NULL) {
  
con <- file(filename,"r")
temp <- readLines(con)
ix <- grep("magicReconstruct",temp)

k <- min(grep("genotype",temp[ix],ignore.case=F))
x <- strsplit(temp[c((ix[k]+2):(ix[k+1]-1))],split=",")
ng <- length(x)
genotypes <- sapply(x,"[",3)

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

write.csv(cbind(map,out1),file=outfile1,row.names=F)
#write.csv(data.frame(id=colnames(out),parent1="F1",parent2="F1"),file="Week5/diaQTL_F2ped.csv",row.names=F)
}