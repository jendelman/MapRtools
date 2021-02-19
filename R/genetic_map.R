#' Multi-point estimation of a genetic map
#' 
#' Multi-point estimation of a genetic map
#' 
#' Uses LOD-score weighted least-squares regression method described by Stam (1993). Markers must be binned (e.g., using \code{\link{LDbin}}) for this function to work properly. Argument \code{n.point} controls how many pairwise distances are used in the linear regression. \code{n.point=2} means only adjacent bins; \code{n.point=3}  means adjacent bins and bins with one intervening marker, etc. Marker names taken from the rownames attribute of \code{x}. 
#' 
#' @param x matrix of pairwise map distances (cM) between the marker-bins for one chromosome
#' @param LOD matrix of LOD scores between marker-bins
#' @param n.point Number of points used for estimation
#' 
#' @return data frame with columns marker,position (in cM)
#' @import CVXR
#' @importFrom Matrix bandSparse
#' @export

genetic_map <- function(x,LOD,n.point=5) {
  
  m <- nrow(x)
  n.point <- min(n.point,m)
  A <- NULL
  y <- NULL
  w <- NULL
  for (j in 2:n.point) {
    n <- m-(j-1)
    A <- rbind(A,bandSparse(n,m,k=c(0,j-1),diagonals=list(rep(-1,n),rep(1,n))))
    y <- c(y,x[cbind(1:n,j:m)])
    w <- c(w,LOD[cbind(1:n,j:m)])
  }
  A <- as.matrix(A)
  stopifnot(all(w!=Inf))
  
  u <- Variable(m)
  objective <- sum(w*(y - A %*% u)^2)
  
  b <- matrix(0,nrow=1,ncol=m)
  b[1] <- 1
  B <- A[1:(m-1),]
  constraints <- list(B%*%u >= 0,b%*%u==0)
  
  prob <- Problem(Minimize(objective),constraints)
  result <- solve(prob)
  ans <- data.frame(marker=rownames(x),position=round(as.vector(result$getValue(u)),2),stringsAsFactors=F)
  return(ans)
}