#' Signatures of selection in S1 populations
#' 
#' Signatures of selection in S1 populations
#' 
#' Genotypes must be coded based on the S1 parental haplotypes, not markers.
#'
#' The null hypothesis is no selection, in which case the expected frequency of genotypes is (AA = 1/4, AB = 1/2, BB = 1/4). Two alternate hypotheses are tested for gametic selection: 1.selection in one sex, 2.selection in both sexes. Two models of zygotic selection are also tested: 1.selection against one homozygote, 2.selection against both homozygotes. The selection coefficient equals the sum of the absolute differences between the observed and expected frequencies. Positive values correspond to selection against A or AA, negative values for selection against B or BB. For zygotic2, positive (negative) values represent selection against (for) homozygotes. 
#' 
#' P-values are computed based on the likelihood ratio test; in other words, the change in deviance is assumed to be chi-squared distributed under the null hypothesis.
#' 
#' @param data data frame with columns: marker, chrom, position, AA, AB, BB. Columns 4-6 have count data.
#' @param alpha significance level
#' 
#' @return list with "plot" and "table" of results:
#' \describe{
#' \item{model}{name of best model}
#' \item{s}{selection coefficient}
#' \item{score}{-log10(p) value}
#' }
#' 
#' @importFrom stats optimize pchisq
#' @import ggplot2
#' @export
#' 

S1_selection <- function(data,alpha=0.05) {
  
  get_x <- function(map) {
    #takes a map with chrom and position and returns x axis values for plotting multiple chromosomes
    a <- tapply(map[,2],map[,1],max)
    n <- length(a)
    m <- tapply(map[,2],map[,1],length)
    b <- c(0,apply(array(1:(n-1)),1,function(k){sum(a[1:k])}))
    x <- map[,2] + rep(b,times=m)
    return(x)
  }
  
  f <- function(s) {
    ifelse(s > 0,s/2,s/6)
  }
  
  freq <- function(s,model) {
    #return zygote frequencies 
    switch(model,
           gametic1=1/4*c(1-s,2,1+s),
           gametic2=1/4*c((1-s)^2,1-s^2,(1+s)^2),
           zygotic1=c(1/4-f(s),1/2+abs(s)/3,1/4-f(-s)),
           zygotic2=1/4*c(1-s,2*(1+s),1-s))
  }
  
  f.LL <- function(s,N,model) {
    p <- freq(s,model)
    return(sum(N*log(p)))
  }
  
  x <- data[,4:6]
  map <- data[,1:3]
  colnames(map) <- c("marker","chrom","position")
  
  models <- c("gametic1","gametic2","zygotic1","zygotic2")
  s <- NULL
  LL <- NULL
  for (i in 1:4) {
    if (i==3) {
      limits <- c(-0.5+1e-3,0.5-1e-3)
    } else {
      limits <- c(-1+1e-3,1-1e-3)
    }
    y <- t(apply(x,1,function(N){
            as.numeric(optimize(f.LL,interval=limits,N=N,model=models[i],maximum=TRUE))}))
    s <- cbind(s,y[,1])
    LL <- cbind(LL,y[,2])
  }
  ix <- apply(LL,1,which.max)
  
  m <- nrow(x)
  LL1 <- LL[cbind(1:m,ix)]
  LL0 <- apply(x,1,function(N){f.LL(s=0,N=N,model="gametic1")})
  score <- -log10(pchisq(q=2*(LL1-LL0),df=1,lower.tail=F))
  output <- data.frame(map,model=models[ix],s=round(s[cbind(1:m,ix)],3),score=round(score,3))
  
  threshold <- -log10(alpha/m)
  plotme <- data.frame(output,color=ifelse(score > threshold,output$model,"NS"))
  chrom <- unique(map$chrom)
  nchr <- length(chrom)
  plotme$model <- factor(plotme$model,levels=c(models,"NS"))
  colours <- c("#66C2A5","#FC8D62","#8DA0CB","#E78AC3","gray60")
  names(colours) <- levels(plotme$model)
  
  if (nchr==1) {
    p <- ggplot(data=plotme,aes(x=.data$position,y=.data$score,colour=.data$color)) + 
      geom_point() + 
      scale_colour_manual(name="Model",values=colours) + theme_bw() + 
      ylab(expression(paste("-log"[10],"(p)"))) + theme(text = element_text(size=13)) + 
      xlab("Position") + geom_hline(yintercept=threshold,linetype=2,col="gray30") 
  } else {
    x <- get_x(map[,2:3])
    break1 <- (tapply(x,map$chrom,max) + tapply(x,map$chrom,min))/2
    break2 <- (tapply(x,map$chrom,max)[-nchr] + tapply(x,map$chrom,min)[-1])/2
    plotme$x <- x
    p <- ggplot(data=plotme,aes(x=.data$x,y=.data$score,colour=.data$color)) +
      geom_point() + 
      scale_colour_manual(name="Model",values=colours) + theme_bw() + 
      scale_x_continuous(name="Chromosome",breaks=break1,labels=chrom,
                         minor_breaks = break2) +
      ylab(expression(paste("-log"[10],"(p)"))) + 
      theme(text = element_text(size=13), 
            panel.grid.major.x = element_blank()) +
      geom_hline(yintercept=threshold,linetype=2,col="gray30")
  }
  return(list(table=output,plot=p))
}
