#' Graphical genotyping
#' 
#' Graphical genotyping
#' 
#' Input matrix \code{geno} should have rownames attribute that matches marker names in the first column of \code{map}. 
#' 
#' @param geno genotype matrix (markers x indiv)
#' @param map data frame with 3 columns (marker, chrom, position)
#' 
#' @return ggplot object 
#' 
#' @export
#' @import ggplot2

plot_haplo <- function(geno,map) {
  
  map$marker <- as.character(map$marker)
  markers <- rownames(geno)
  stopifnot(markers %in% map$marker)
  map <- map[map$marker %in% markers,]
  n <- ncol(geno)
  
  chroms <- unique(map$chrom)
  n.chrom <- length(chroms)
  plot.data <- NULL
  for (i in 1:n.chrom) {
    ix <- which(map$chrom==chroms[i])
    m <- length(ix)
    x <- map$position[ix]
    d <- diff(x)/2
    xmin <- c(x[1]-d[1],x[-m]+d)
    xmax <- c(x[-m]+d,x[m]+d[m-1])
    plot.data <- rbind(plot.data,data.frame(z=factor(as.vector(geno[ix,])),chrom=chroms[i],xmin=xmin,xmax=xmax,ymin=rep(1:n-1,each=m),ymax=rep(1:n,each=m)))
  }
  plot.data$chrom <- factor(plot.data$chrom)
  p <- ggplot(data=plot.data) + geom_rect(aes(fill=z,xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax)) + theme_bw() + scale_fill_manual(name="Dosage",values=c("red","green","blue"),na.value="grey50") + xlab("Position") + ylab("Sample Number")
  
  if (n.chrom > 1) {
    p <- p + facet_wrap(~chrom,scales = "free_x")
  }
  return(p)
}
