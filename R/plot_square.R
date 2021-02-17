#' Plot square (dis)similarity matrix
#' 
#' Plot square (dis)similarity matrix
#'
#' Can be used to plot squared correlation, recomb frequency, LOD and more. By default, \code{lims} equals (0,median,max) 
#'
#' @param data squared correlation matrix
#' @param lims numeric 3-vector with the low,mid,high points for the colors
#' 
#' @return ggplot2 variable
#' @export
#' @import ggplot2

plot_square <- function(data,lims=NULL) {
  z <- as.vector(data)
  z <- ifelse(z==Inf,NA,z)
  if (is.null(lims)) {
    lims <- c(0,mean(z,na.rm=T),max(z,na.rm=T))
  }
  stopifnot(length(lims)==3)
  m <- nrow(data)
  data2 <- expand.grid(x=1:m,y=1:m)
  data2$z <- z  #vectorize the matrix
  ggplot(data2,aes(x=x,y=y)) + geom_tile(aes(fill=z)) + xlab("") + ylab("") +
    scale_fill_gradient2(low="yellow",mid="orange",midpoint=lims[2],high="red",name="",limits=lims[c(1,3)]) + coord_fixed(ratio=1)
  
}
