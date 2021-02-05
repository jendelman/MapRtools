#' Visualize haplotype dosage
#' 
#' Visualize haplotype dosage in diploid biparental population from two inbreds
#' 
#' Input matrix \code{geno} should have rownames attribute that matches marker names in the first column of \code{map}. 
#' 
#' @param geno matrix of haplotype dosages (markers x indiv)
#' @param map data frame with 3 columns (marker, chrom, position)
#' 
#' @return ggplot object 
#' 
#' @export
#' @import ggplot2

plot_haplo <- function(geno,map) {
  
  map$marker <- as.character(map$marker)
  map <- map[map$marker %in% rownames(geno),]
  geno <- geno[match(map$marker,rownames(geno)),]
  n <- ncol(geno)
  m <- nrow(geno)
  x <- map$position
  xmin <- c(0,x[-m] + diff(x)/2)
  xmax <- c(xmin[-1],x[m])
  
  plotme <- data.frame(z=factor(as.vector(geno)),chrom=factor(map$chrom),xmin=xmin,xmax=xmax,ymin=rep(1:n-1,each=m),ymax=rep(1:n,each=m))
  p <- ggplot(data=plotme) + geom_rect(aes(fill=z,xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax)) + theme_bw() + scale_fill_manual(name="Dosage",values=c("red","green","blue"),na.value="grey50") + xlab("Position") + ylab("Sample Number")
  n.chrom <- length(levels(plotme$chrom))
  if (n.chrom > 1) {
    p <- p + facet_wrap(~chrom,scales = "free_x")
  }
  return(p)
}
