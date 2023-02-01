#' Plot and filter markers based on genotype frequency vs position
#' 
#' Plot and filter markers based on genotype frequency vs position
#' 
#' Genotypes should be coded 0,1,2. Markers are removed if their residual to the fitted spline exceeds \code{thresh}. Markers are assumed to be ordered. Function designed to be used for one chromosome.
#' 
#' @param geno haplotype dosage matrix (markers x indiv)
#' @param thresh threshold for removing markers (see Details)
#' @param span parameter to control degree of smoothing for spline (higher = less smooth)
#' 
#' @return List containing
#' \describe{
#' \item{outliers}{character vector of marker names}
#' \item{plot}{ggplot2 variable}
#' }
#' @export
#' @importFrom stats loess
#' @import ggplot2

plot_genofreq <- function(geno,thresh=0.1,span=0.3) {
  m <- nrow(geno)
  n <- ncol(geno)
  gf <- apply(geno,1,function(x){table(factor(x,levels=0:2))})/n
  outliers <- integer(0)
  lans <- vector("list",3)
  for (i in 1:3) {
    lans[[i]] <- loess(y~x,data.frame(y=gf[i,],x=1:m),span=span)
    outliers <- union(outliers,which(abs(lans[[i]]$residuals) > thresh))
  }
  keep <- setdiff(1:m,outliers)
  
  plotme <- data.frame(y1=c(gf[1,],gf[2,],gf[3,]),y2=c(lans[[1]]$fitted,lans[[2]]$fitted,lans[[3]]$fitted),x=rep(1:m,3),genotype=factor(rep(c("0","1","2"),each=m)),outlier=factor(ifelse(1:m %in% outliers,1,2)))
  col <- c("red","black")
  names(col) <- 1:2
  p <- ggplot(data=plotme,aes(x=x,colour=outlier)) + geom_point(mapping=aes(y=y1)) + facet_wrap(~genotype,ncol=1,nrow=3) + theme_bw() + xlab("") + ylab("Frequency") + geom_line(aes(y=y2),col="blue") + scale_color_manual(name="",values=col) + guides(colour="none")
  return(list(outliers=rownames(geno)[outliers],plot=p))
}