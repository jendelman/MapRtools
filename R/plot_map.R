#' Plot data against position and curate based on trend
#' 
#' Plot data against position and curate based on trend
#' 
#' For data from a single chromosome, markers are removed when their residual to the fitted spline exceeds \code{thresh}. If NULL, no curation occurs.
#' 
#' @param data data frame with 3 columns: chrom, position, plotting variable
#' @param thresh threshold for removing markers (see Details)
#' @param span smoothing parameter for loess spline (higher = less smooth)
#'
#' @return ggplot
#'
#' @return List containing
#' \describe{
#' \item{outliers}{outlier positions}
#' \item{plot}{ggplot2 variable}
#' }
#' @export
#' @importFrom stats loess
#' @import ggplot2

plot_map <- function(data,thresh=NULL,span=0.3) {
  
  get_x <- function(map) {
    #takes a map with chrom and position and returns x axis values for plotting multiple chromosomes
    a <- tapply(map[,2],map[,1],max)
    n <- length(a)
    m <- tapply(map[,2],map[,1],length)
    b <- c(0,apply(array(1:(n-1)),1,function(k){sum(a[1:k])}))
    x <- map[,2] + rep(b,times=m)
    return(x)
  }
  trait <- colnames(data)[3]
  m <- nrow(data)
  colnames(data) <- c("chrom","pos","y")
  data$chrom <- as.character(data$chrom)
  chrom <- unique(data$chrom)
  data$pos <- as.integer(data$pos)
  nchr <- length(unique(data$chrom))
  data$color <- factor(ifelse(as.integer(factor(data$chrom))%%2==1,1,0))

  outliers <- integer(0)
  
  if (nchr==1L) {
    lans <- loess(y~x,data.frame(y=data$y,x=data$pos),span=span)
    data$fitted <- lans$fitted
    
    if (!is.null(thresh)) {
      data$outlier <- ifelse(abs(lans$residuals) > thresh,"1","0")
      outliers <- which(data$outlier=="1")
      col <- c("0"="black","1"="red")
      p <- ggplot(data=data,aes(x=pos,colour=outlier)) + 
        geom_point(mapping=aes(y=y)) + 
        geom_line(mapping=aes(y=fitted),col="blue",linetype=2) + 
        theme_bw() + xlab("") + ylab(trait) + 
        scale_color_manual(name="",values=col) + guides(colour="none")
    } else {
      p <- ggplot(data=data,aes(x=pos)) + 
        geom_point(mapping=aes(y=y)) + 
        geom_line(mapping=aes(y=fitted),col="blue",linetype=2) + 
        theme_bw() + xlab("") + ylab(trait)
    }
    
  } else {
    data$x <- get_x(data[,1:2])
    break1 <- (tapply(data$x,data$chrom,max) + tapply(data$x,data$chrom,min))/2
    break2 <- (tapply(data$x,data$chrom,max)[-nchr] + tapply(data$x,data$chrom,min)[-1])/2
    p <- ggplot(data=data,aes(x=.data$x,y=.data$y,colour=.data$color)) +
      geom_point() + 
      scale_colour_manual(values=c("#21908c","#440154")) + theme_bw() + 
      scale_x_continuous(name="Chromosome",breaks=break1,labels=chrom,
                         minor_breaks = break2) +
      theme(text = element_text(size=13), 
            panel.grid.major.x = element_blank()) + guides(colour="none") + ylab("")
  }
  return(list(outliers=outliers,p=p))
}
