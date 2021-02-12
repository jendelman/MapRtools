#' Trim a linkage group based on genotype frequencies
#' 
#' Trim a linkage group based on genotype frequencies
#' 
#' This function should only be run on a single linkage group (to form the linkage groups, use \code{\link{LG}}. If \code{thresh} is a numeric vector with multiple LOD thresholds, the function returns a plot showing the impact of the threshold on genotype frequencies. If \code{thresh} is a single value, the function returns a vector of the marker names that are retained. The rownames of \code{geno} and \code{LODmat} must match.
#' 
#' @param geno matrix of haplotype dosages (markers x samples)
#' @param LODmat matrix of LOD scores for the markers
#' @param thresh numeric vector of thresholds for clusterings
#' 
#' @return Either a ggplot2 object or a vector of marker names (see Details)
#' @export
#' @import ggplot2
#' @importFrom stats hclust cutree

LGtrim <- function(geno,LODmat,thresh) {
  
  dLOD <- exp(-LODmat)
  diag(dLOD) <- 0
  x <- hclust(as.dist(dLOD),method="single")
  thresholds <- exp(-thresh)
  tmp <- cutree(x,h=min(thresholds))
  max.clust <- max(tmp)
  n.thresh <- length(thresholds)
  data <- NULL
  for (i in 1:n.thresh) {
    tmp <- cutree(x,h=thresholds[i])
    tab <- table(tmp)
    ix <- tmp==as.integer(names(tab)[which.max(tab)])
    geno.freq <- apply(geno[ix,],1,function(x){
      tab <- table(factor(x,levels=0:2))
      tab/sum(tab)
    })
    data <- rbind(data,data.frame(y=as.vector(geno.freq),x=factor(rep(0:2,ncol(geno.freq))),thresh=thresh[i]))
  }
  if (n.thresh==1) {
    return(rownames(geno)[ix])
  } else {
    data$thresh <- factor(data$thresh)
    p <- ggplot(data,aes(x=x,y=y)) + geom_violin() + facet_wrap(~thresh) + theme_bw() + ylab("Frequency") + xlab("Dosage")
    return(p)
  }
}
