#' Make linkage groups based on clustering
#' 
#' Make linkage groups based on clustering
#' 
#' If \code{thresh} is a numeric vector with multiple LOD thresholds, the function returns a plot showing the number of markers per LG. If \code{thresh} is a single value, the function returns a data frame with the LG assignment for each marker. LGs are numbered from the largest to smallest group. 
#' 
#' @param LODmat matrix of LOD scores for the marker bins
#' @param thresh numeric vector of thresholds for clusterings
#' 
#' @return Either a ggplot2 object or data frame of linkage groups (see Details)
#' @export
#' @import ggplot2
#' @importFrom stats hclust cutree

LG <- function(LODmat,thresh=seq(2,20,by=2)) {
  
  dLOD <- exp(-LODmat)
  diag(dLOD) <- 0
  x <- hclust(as.dist(dLOD),method="single")
  thresholds <- exp(-thresh)
  tmp <- cutree(x,h=min(thresholds))
  max.clust <- max(tmp)
  n.thresh <- length(thresholds)
  LGans <- matrix(0,nrow=n.thresh,ncol=max.clust)
  rownames(LGans) <- thresholds
  for (i in 1:n.thresh) {
    tmp <- cutree(x,h=thresholds[i])
    LGans[i,1:max(tmp)] <- sort(table(tmp),decreasing=T)
  }
  if (n.thresh==1) {
    out <- data.frame(marker=rownames(LODmat),LG=as.integer(tmp),stringsAsFactors = F)
    tab <- table(out$LG)
    nG <- length(tab)
    tmp <- match(1:nG,order(tab,decreasing=T))
    out$LG <- tmp[out$LG]
    return(out)
  } else {
    rownames(LGans) <- thresh
    plot.data <- data.frame(count=as.vector(t(LGans)),thresh=factor(rep(thresh,each=max.clust)),
                          stripe=factor(rep(1:max.clust,times=n.thresh),ordered=T,levels=seq(max.clust,1,-1)))
  
    stripe.color <- rep(c("red","white"),ceiling(max.clust/2))[1:max.clust]
    p <- ggplot(plot.data,aes(y=count,x=thresh,fill=stripe)) + geom_col(position="stack",colour="black") + scale_fill_manual(values=stripe.color) + guides(fill="none") + theme_bw() + xlab("Threshold")
    return(p)
  }
}
