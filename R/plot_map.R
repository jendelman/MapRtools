#' Plots data against map
#' 
#' Plots data against map
#' 
#' @param data data frame with 3 columns: chrom, position, y (the plotting variable)
#' 
#' @return ggplot
#'
#' @import ggplot2
#' @export
#' 

plot_map <- function(data) {
  
  get_x <- function(map) {
    #takes a map with chrom and position and returns x axis values for plotting multiple chromosomes
    a <- tapply(map[,2],map[,1],max)
    n <- length(a)
    m <- tapply(map[,2],map[,1],length)
    b <- c(0,apply(array(1:(n-1)),1,function(k){sum(a[1:k])}))
    x <- map[,2] + rep(b,times=m)
    return(x)
  }
  colnames(data) <- c("chrom","pos","y")
  data$chrom <- as.character(data$chrom)
  chrom <- unique(data$chrom)
  data$pos <- as.integer(data$pos)
  nchr <- length(unique(data$chrom))
  data$color <- factor(ifelse(as.integer(factor(data$chrom))%%2==1,1,0))
  
  if (nchr==1) {
    data$x <- data$pos/1e6 
    p <- ggplot(data=data,aes(x=.data$x,y=.data$y)) + 
      geom_point() + theme_bw() +
      theme(text = element_text(size=13)) + 
      xlab("Position (Mb)") 
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
  return(p)
}
