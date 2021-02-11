#' Create marker bins based on LD
#' 
#' Create marker bins based on LD
#' 
#' Bins are created based on hierarchical clustering with \code{hclust} and \code{method='single'}, using \eqn{1-r^2} as the dissimilarity metric. The argument \code{r2.thresh} controls the height for cutting the dendrogram to create the bins. The marker with the least missing data for each bin is chosen to represent it.
#' 
#' @param geno matrix of haplotype dosages (markers x indiv)
#' @param r2.thresh threshold for binning
#' 
#' @return List containing
#' \describe{
#' \item{bins}{data frame with two columns: marker,bin}
#' \item{geno}{genotype matrix for the bins}
#' \item{r2}{r2 matrix for the bins}
#' }
#' 
#' @export
#' @importFrom stats hclust cor as.dist cutree

LDbin <- function(geno,r2.thresh=0.99) {
  markers <- rownames(geno)
  m <- length(markers)
  stopifnot(m>1)
    
  r2 <- cor(t(geno),use = "pairwise.complete")^2
  clustans <- hclust(as.dist(1-r2),method="single")
  bins <- cutree(clustans,h=1-r2.thresh) #vector of bin numbers
  bins <- data.frame(marker=names(bins),bin=as.integer(bins))
  x <- split(bins$marker,bins$bin)
  nmiss <- apply(geno,1,function(z){sum(is.na(z))})
  y <- split(nmiss,bins$bin)
  bin.marker <- mapply(FUN=function(x,y){x[which.min(y)]},x=x,y=y)
  return(list(bins=bins,geno=geno[bin.marker,],r2=r2[bin.marker,bin.marker]))
}  
