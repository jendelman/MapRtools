setwd("~/PanGenome/DMv6.1/")
x <- read.table("DMchr_end.txt",as.is=T)[,1]
start <- c(1,x[1:11]+1)
bp.length <- x-start+1 
con <- file("DM_1-3_516_R44_potato_genome_assembly.v6.1.fa","r")
chroms <- c("chr01","chr02","chr03","chr04","chr05","chr06","chr07","chr08","chr09","chr10","chr11","chr12")
for (i in 1:12) {
  print(i)
  con2 <- file(paste(chroms[i],"fa",sep="."),"w")
  temp <- readLines(con,n=bp.length[i])
  writeLines(text=temp,con=con2)
  close(con2)
}
close(con)

get_context <- function(chrom,bp,flank=100) {
  con <- file(paste(chrom,"fa",sep="."),"r")
  start.c <- (bp-flank) %% 60
  if (start.c==0) {
    start.c <- 60
    start.L <- (bp-flank) %/% 60  
  } else {
    start.L <- (bp-flank) %/% 60 + 1
  }
  
  end.c <- (bp+flank) %% 60
  if (end.c==0) {
    end.c <- 60
    end.L <- (bp+flank) %/% 60  
  } else {
    end.L <- (bp+flank) %/% 60 + 1
  }
  deltaL <- end.L - start.L
  tmp <- readLines(con,n=start.L) #includes header
  tmp <- readLines(con,n=1) 
  if (deltaL==0) {
    out <- substr(tmp,start.c,end.c)
    return(out)
  } else {
    out <- substr(tmp,start.c,60)
  }
  if (deltaL > 1) {
    tmp <- readLines(con,n=deltaL-1)
    out <- paste(c(out,tmp),collapse="")
  }
  tmp <- readLines(con,n=1)
  out <- paste(out,substr(tmp,1,end.c),sep="")
  close(con)
  return(out)
}

reverse_complement <- function(x) {
  n <- nchar(x)
  out <- character(n)
  for (i in 1:n) {
    y <- substr(x,i,i)
    k <- match(y,c("A","C","G","T"))
    out[n-i+1] <- switch(k,"T","G","C","A")
  }
  return(paste(out,collapse=""))
}
reverse_complement("ACGTTTC")

data <- read.csv("Round2_design.csv",as.is=T)
m <- nrow(data)
checkit <- logical(m)
out <- character(m)
for (i in 1:m) {
  print(i)
  out[i] <- get_context(chrom=data$Chrom[i],bp=data$Position[i],flank=150)
  tmp <- data$TargetSequence[i]
  john <- paste(c(substr(tmp,1,50),substr(tmp,52,52),substr(tmp,56,105)),collapse="")
  checkit[i] <- john==substr(out[i],101,201)
}
data$LongTarget <- out

Ryadg <- "ACGTGCTAACTAGTTAGGGATTCAAATTCAMRATTGTATTAAACCCGGATATACATATACAGGGAAGYTTTAACCACACATGMAAGGTTCAGATATCCAWG"
nchar(Ryadg)
tmp <- get_context(chrom="chr11",bp=2499552-50,flank=150)
snpST00073 <- paste(c(substr(tmp,1,150),"[T/C]",substr(tmp,152,301)),collapse="")

Rysto <- "CTCAAGCGGAATAACCCCTTTGGCTTCTAGAGTTGAACGTTTC"
nchar(Rysto)
tmp <- get_context(chrom="chr12",bp=2352763-21,flank=150)
substr(tmp,151,151)
reverse_complement(tmp)==Rysto
snpST00082 <- paste(c(substr(tmp,1,150),"[A/C]",substr(tmp,152,301)),collapse="")
snpST00082

kasp <- data.frame(snpID=c("snpST00073","snpST00082"),TargetSequence=c(snpST00073,snpST00082),Chrom=c("chr11"
,"chr12"),Position=c(2499552-50,2352763-21),Alleles=c("[T/C]","[A/C]"))
write.csv(kasp,"Ry_kasp.csv")
