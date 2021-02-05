#' Plot heatmap of LD
#' 
#' Plot heatmap of LD (squared correlation)
#'
#' @param r2 squared correlation matrix
#' 
#' @return ggplot2 variable
#' @export

plot_r2 <- function(r2) {
  m <- nrow(r2)
  data <- expand.grid(x=1:m,y=1:m)
  data$z <- as.vector(r2)  #vectorize the matrix
  ggplot(data,aes(x=x,y=y)) + geom_tile(aes(fill=z)) + xlab("") + ylab("") +
    scale_fill_gradient2(low="yellow",mid="orange",midpoint=0.5,high="red",name=expression(r^2),limits=c(0,1)) + coord_fixed(ratio=1)
}
