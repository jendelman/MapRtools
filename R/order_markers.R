#' Order markers by solving the TSP
#' 
#' Order markers by solving the TSP
#' 
#' Uses R package seriation to minimize the distance between adjacent markers. For example, \code{x} could be a matrix of recombination frequencies or monotone decreasing transformation of LOD scores.
#'
#' @param x distance matrix
#' 
#' @return a list containing
#' \describe{
#' \item{path}{optimized order as a vector of integers}
#' \item{distance}{sum of adjacent distances}
#' }
#' 
#' @export
#' @importFrom seriation seriate get_order

order_markers <- function(x) {
  o <- seriate(as.dist(x),method="TSP")
  path <- get_order(o)
  x_ordered <- x[path,path]
  m <- nrow(x)
  d <- x_ordered[cbind(1:(m-1),2:m)]
  return(list(path=path,distance=sum(d)))
}
