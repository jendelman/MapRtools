#' Plot marker coverage of the genome
#' 
#' Plot marker coverage of the genome
#' 
#' Distances can be given in Mb or cM
#' 
#' @param map Data frame with two columns: chrom, position 
#' @param limits Data frame with two columns: chrom, length
#' 
#' @return ggplot2 variable
#' 
#' @import ggplot2
#'
plot_coverage <- function(map,limits) {
  colnames(map) <- c("chrom","position")
  chroms <- unique(map$chrom)
  n.chrom <- length(chroms)
  k <- match(map$chrom,chroms)
  y <- c(k-0.1,1:n.chrom)
  yend <- c(k+0.1,1:n.chrom)
  x <- c(map$position,rep(0,n.chrom))
  xend <- c(map$position,limits$length)
  p <- ggplot(data=data.frame(x=x,y=y,xend=xend,yend=yend),aes(x=x,y=y,xend=xend,yend=yend)) +  geom_segment() + theme_bw() + xlab("Position") + scale_y_continuous(name="Chromosome",breaks=1:n.chrom,labels=chroms,minor_breaks=NULL)
  return(p)
}
