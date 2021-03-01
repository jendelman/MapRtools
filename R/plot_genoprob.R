#' Plot genotype probabilities for one chromosome
#' 
#' Plot genotype probabilities for one chromosome
#' 
#' Names for the genotypes are taken from the colnames of \code{genoprob}.
#' 
#' @param genoprob matrix (markers x genotypes) of probabilities for one individual
#' @param map map data frame (markers,chrom,position)
#' 
#' @return ggplot object 
#' 
#' @export
#' @import ggplot2

plot_genoprob <- function(genoprob,map) {

  colnames(map) <- c("marker","chrom","position")
  map$marker <- as.character(map$marker)
  markers <- rownames(genoprob)
  stopifnot(markers %in% map$marker)
  
  map <- map[map$marker %in% markers,]
  x <- map$position
  m <- length(x)
  xmin <- c(0,x[-m] + diff(x)/2)
  xmax <- c(xmin[-1],x[m])
  
  ng <- ncol(genoprob)
  y1 <- 0:(ng-1)
  
  plotme <- data.frame(z=as.vector(genoprob),xmin=xmin,xmax=xmax,ymin=rep(y1,each=m),ymax=rep(y1+1,each=m))
  p <- ggplot(data=plotme) + geom_rect(aes(fill=z,xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax)) + theme_bw() + scale_fill_distiller(name="Probability",palette="Blues",direction=1) + scale_y_continuous(name="Genotype",labels=colnames(genoprob),breaks=(1:ng)-0.5) + xlab("Position")
  
  return(p)
}
