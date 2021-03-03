#' Make RABBIT input files for diploid diallel population
#'
#' Make RABBIT input files for diploid diallel population
#' 
#' Populations must be numbered in \code{ped} corresponding to their position in \code{geno}. Founders are not included in \code{ped}. All genotype matrices must have identical markers. Genetic map position should be in cM. Genotypes need to be coded according to RABBIT format. 
#'
#' @param ped data frame with pedigree (pop,parent1,parent2)
#' @param geno list of genotype matrices (markers x indiv), one for each population in \code{ped}
#' @param geno.founder matrix of genotype data for the founders (markers x indiv)
#' @param map genetic map (marker,chromosome,position)
#' @param outstem name for output files
#' 
#' @export
#' @importFrom utils write.table

rabbit_diallel <- function(ped,geno,geno.founder,map,outstem) {
  colnames(ped) <- c("id","parent1","parent2")
  ped$parent1 <- as.character(ped$parent1)
  ped$parent2 <- as.character(ped$parent2)
  founders <- sort(union(ped$parent1,ped$parent2))
  founders <- setdiff(founders,"0")
  stopifnot(founders %in% colnames(geno.founder))
  
  nf <- length(founders)
  if (nf > 1) {
    geno.founder <- geno.founder[,founders]
    colnames(geno.founder) <- founders
  }
  funnel <- paste(1:nf,collapse="-")
  
  npop <- length(geno)
  ped2 <- data.frame(Generation=rep(0,nf),MemberID=1:nf,Gender=rep(0,nf),MotherID=rep(0,nf),FatherID=rep(0,nf))
  ped2 <- rbind(ped2,data.frame(Generation=rep(1,npop),MemberID=1:npop+nf,Gender=rep(0,npop),MotherID=match(ped$parent1,founders),FatherID=match(ped$parent2,founders)))
  
  id <- lapply(geno,colnames)
  
  ped3 <- NULL
  geno2 <- t(geno.founder)
  for (i in 1:npop) {
    idi <- id[[i]]
    n <- length(idi)
    tmp <- data.frame(OffspringID=idi,MemberID=rep(i+nf,n),Funnelcode=rep(funnel,n))
    ped3 <- rbind(ped3,tmp)
    geno2 <- rbind(geno2,t(geno[[i]]))
  }
  marker <- colnames(geno2)
  stopifnot(all(marker %in% map$marker))
  
  m <- length(marker)
  map <- map[match(marker,map$marker),]
  geno3 <- rbind(nfounder=c(nf,rep("",m-1)),t(map),geno2)
  geno3[4,] <- gsub(" ","",geno3[4,]) #removes whitespace
  write.table(x=geno3,file=paste(outstem,"rabbit_geno.csv",sep=""),quote = FALSE,na="NN",col.names=FALSE,sep=",")
  
  ped4a <- rbind(c("Pedigree-Information","DesignPedigree",rep("",3)),colnames(ped2),ped2)
  colnames(ped4a) <- NULL
  ped4b <- rbind(c("Pedigree-Information","SampleInfo",rep("",3)),colnames(ped3),cbind(ped3,rep("",nrow(ped3)),rep("",nrow(ped3))))
  colnames(ped4b) <- NULL
  ped4 <- rbind(as.matrix(ped4a),as.matrix(ped4b))
  write.table(x=ped4,file=paste(outstem,"rabbit_ped.csv",sep=""),quote = FALSE,na="NN",col.names=FALSE,sep=",",row.names=F)
}