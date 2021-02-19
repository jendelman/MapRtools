#' Multi-point estimation of a genetic map
#' 
#' Multi-point estimation of a genetic map
#' 
#' Argument \code{n.point} controls how many pairwise distances are used in the linear regression. \code{n.point=2} means only adjacent bins; \code{n.point=3}  means adjacent bins and bins with one intervening marker, etc. Marker names taken from the rownames attribute of \code{x}. For multi-point estimation to be effective, \code{x} should be based on marker-bins. 
#' 
#' @param x matrix of pairwise map distances (cM) between the marker-bins for one chromosome
#' @param n.point Number of points used for estimation
#' 
#' @return data frame with columns marker,position (in cM)
#' @import CVXR
#' @importFrom Matrix bandSparse
#' @export

genetic_map <- function(x,n.point=3) {
  
  m <- nrow(x)
  n.point <- min(n.point,m)
  A <- NULL
  a <- NULL
  for (j in 2:n.point) {
    n <- m-(j-1)
    A <- rbind(A,bandSparse(n,m,k=c(0,j-1),diagonals=list(rep(-1,n),rep(1,n))))
    a <- c(a,x[cbind(1:n,j:m)])
  }
  A <- as.matrix(A)
  
  y <- Variable(m)
  objective <- sum((a - A %*% y)^2)
  
  b <- matrix(0,nrow=1,ncol=m)
  b[1] <- 1
  B <- A[1:(m-1),]
  constraints <- list(B%*%y >= 0,b%*%y==0)
  
  prob <- Problem(Minimize(objective),constraints)
  result <- solve(prob)
  ans <- data.frame(marker=rownames(x),position=round(as.vector(result$getValue(y)),2),stringsAsFactors=F)
  return(ans)
}