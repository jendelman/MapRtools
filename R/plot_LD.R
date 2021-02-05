#' Plot LD vs distance
#' 
#' Plot LD vs distance
#' 
#' A monotone decreasing, convex spline is fit using R package \code{scam}. The input matrix \code{r2} should have rownames attribute that matches marker names in the first column of \code{map}. 
#' 
#' @param r2 squared correlation matrix
#' @param map data frame with 3 columns (marker, chrom, position)
#' @param max.pair maximum number of r2 pairs for the spline
#' @param dof degrees of freedom for the spline
#' 
#' @return ggplot object 
#' 
#' @export
#' @import ggplot2
#' @import scam
#' @importFrom stats dist

plot_LD <- function(r2,map,max.pair=1e4,dof=5) {
  
  map$marker <- as.character(map$marker)
  map <- map[map$marker %in% rownames(r2),]
  ix <- match(map$marker,rownames(r2))
  r2 <- r2[ix,ix]
  
  chroms <- unique(map$chrom)
  n.chrom <- length(chroms)
  result <- NULL
  for (i in 1:n.chrom) {
    ix <- which(map$chrom==chroms[i])
    m <- length(ix)
    tmp <- expand.grid(col=1:m,row=1:m)
    tmp <- tmp[tmp$row >= tmp$col,]  #only need lower triangular
    r2.vec <- as.vector(r2[ix,ix][cbind(tmp$row,tmp$col)])
    
    d <- as.matrix(dist(matrix(map$position[ix],ncol=1))) #distance matrix
    d.vec <- as.vector(d[cbind(tmp$row,tmp$col)])
    
    result <- rbind(result,data.frame(d=d.vec,r2=r2.vec))
  }
  
  #Hexbin plot
  p <- ggplot(data=result,aes(x=d,y=r2)) + stat_binhex(mapping=aes(colour=..count..)) + ylab(expression(r^2)) + xlab("Distance") + theme_bw() + ylim(0,1)
  
  #Spline
  scam.ans <- scam(formula=r2~s(d,bs=c("mdcx"),k=dof),data=result[sample(1:nrow(result),max.pair),])
  dmax <- max(result$d)
  predans <- predict.scam(scam.ans,newdata=data.frame(d=seq(0,dmax,length.out = 500)))
  p <- p + geom_line(data=data.frame(d=seq(0,dmax,length.out = 500),r2=predans),colour="red")
  return(p)
}
