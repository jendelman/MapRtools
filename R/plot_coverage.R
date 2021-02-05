#' Plot marker coverage of the genome
#' 
#' Plot marker coverage of the genome
#' 
#' If \code{limits} not provided, then the maximum values in \code{map} are used.
#' 
#' @param map data frame with columns chrom & position
#' @param limits optional data frame with columns chrom & position, with the maximum length for each chromosome
#' 
#' @return ggplot2 variable
#' 
#' @export
#' @import ggplot2
#'
plot_coverage <- function(map,limits=NULL) {
  if (is.null(limits)) {
    tmp <- tapply(map$position,map$chrom,max)
    limits <- data.frame(chrom=names(tmp),position=as.numeric(tmp))
  }
  map$chrom <- as.character(map$chrom)
  chroms <- unique(map$chrom)
  n.chrom <- length(chroms)
  k <- match(map$chrom,chroms)
  y <- c(k-0.1,1:n.chrom)
  yend <- c(k+0.1,1:n.chrom)
  x <- c(map$position,rep(0,n.chrom))
  xend <- c(map$position,limits$position)
  p <- ggplot(data=data.frame(x=x,y=y,xend=xend,yend=yend),aes(x=x,y=y,xend=xend,yend=yend)) +  geom_segment() + theme_bw() + xlab("Position") + scale_y_continuous(name="Chromosome",breaks=1:n.chrom,labels=chroms,minor_breaks=NULL)
  return(p)
}
