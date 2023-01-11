## Constuct a table with coding-noncoding, noncoding-noncoding, coding-anticoding
# Id1,id2, promoter distance, overlap length, single1, single2, single nuclear, single cell3 
# load sense-antisense gene pairs extracted from araport11
co_non <- read.table("/Users/xung0001/Desktop/coding-noncoding.txt",head=T,sep = ",")
co_co <- read.table("/Users/xung0001/Desktop/coding_anticoding.txt",head=T,sep = ",")
non_non <- read.table("/Users/xung0001/Desktop/noncoding_noncoding.txt",head=T,sep = ",")
# load single cell tables
sin1 <- read.csv("/Users/xung0001/Downloads/GSE46226_SC_Expression.csv")
row.names(sin1) <- sin1$Locus
sin1 <- sin1[,-1]
sin1 <- sin1[,-(4:7)]
head(sin1)
dim(sin1)
sin1<- sin1[rowSums(sin1[,1:27])>0,]
sin2 <- read.csv("/Users/xung0001/Downloads/GSE74488_sc_expression.csv")
row.names(sin2) <- sin2$Locus
sin2 <- sin2[,-1]
sin2 <- sin2[,1:30]
dim(sin2)
sin2<- sin2[rowSums(sin2[,1:30])>0,]
head(sin2)
sin3 <- read.csv("/Users/xung0001/Downloads/GSE155304_sCellRNA-seq_avg_UMI_counts_per_cluster.tsv",sep="\t")
row.names(sin3) <- sin3$gene
head(sin3)
dim(sin3)
sin3<- sin3[rowSums(sin3[,2:21])>0,]
snuc <- read.csv("/Users/xung0001/Downloads/GSE155304_sNucRNA-seq_avg_UMI_counts_per_cluster.tsv",sep="\t")
row.names(snuc) <-snuc$gene
head(snuc)
dim(snuc)
snuc<- snuc[rowSums(snuc[,2:21])>0,]
# load genome information from sense file 
library(GenomicRanges)
# open adjusted araport11 Large Granges from desktop from sebastian group 
genes_araport_adj <- readRDS("genes_araport_adj.RDS")
unique(factor(genes_araport_adj$tx_type))
##load the source of coding-noncoding corresponding sense mRNA and antisense noncoding
# extract subgroup Granges and their IDs from  genes_araport_adj
senslncRNA <- genes_araport_adj[which(genes_araport_adj$tx_type=="lnc_RNA"),]
antilncRNA <- genes_araport_adj[which(genes_araport_adj$tx_type=="antisense_lncRNA"),]
lnc_rna <- c(senslncRNA,antilncRNA)
mRNA <- genes_araport_adj[which(genes_araport_adj$tx_type=="mRNA"),]
tRNA <- genes_araport_adj[which(genes_araport_adj$tx_type=="tRNA"),]
ncRNA <- genes_araport_adj[which(genes_araport_adj$tx_type=="ncRNA"),]
rRNA <- genes_araport_adj[which(genes_araport_adj$tx_type=="rRNA"),]
snoRNA <- genes_araport_adj[which(genes_araport_adj$tx_type=="snoRNA"),]
snRNA <- genes_araport_adj[which(genes_araport_adj$tx_type=="snRNA"),]
noncodingRNA <- c(tRNA,ncRNA,rRNA,snoRNA,snRNA,lnc_rna)
noncodingRNA1 <- c(tRNA,ncRNA,rRNA,snoRNA,snRNA,antilncRNA)
noncodingRNA2 <- c(tRNA,ncRNA,rRNA,snoRNA,snRNA)
noncodingRNA3 <- c(ncRNA,rRNA,snoRNA,snRNA)
noncodingRNA4 <- c(rRNA,snoRNA,snRNA)
noncodingRNA5 <- c(rRNA,snRNA)
# generate dataframe contains mRNA and noncodingRNA genome information for coding-noncoding 
mRNA <- as.data.frame(mRNA)
mRNA <- mRNA[,1:6]
noncodingRNA <- as.data.frame(noncodingRNA)
noncodingRNA <- noncodingRNA[,1:6]
s <- rbind(mRNA,noncodingRNA)
dim(s[which(co_non$id1%in% s$gene_id),])
dim(s[which(co_non$id2%in% s$gene_id),])

dim(mRNA[which(co_co$V1%in% mRNA$gene_id),])
dim(mRNA[which(co_co$V2%in% mRNA$gene_id),])
## add promoter distance 
non <- co_non
colnames(non)[1:2] <- c("id1","id2")
head(non)
row.names(s) <- s$gene_id
length <- c()
for (i in 1:nrow(non))
{print(i)
  if(s[non[i,2],]$strand=="-")
  {length <-rbind(length,data.frame(id1=non[i,1],id2=non[i,2], length=s[non[i,2],]$end-s[non[i,1],]$start))}
  if(s[non[i,2],]$strand=="+") {length <-rbind(length,data.frame(id1=non[i,1],id2=non[i,2], length=s[non[i,1],]$end-s[non[i,2],]$start))
  }
  
}
non$promoterdis <- length$length
non <- non[order(non$promoterdis,decreasing = F),]
non <- non[-which(non$promoterdis==0),]
## add gene length
length1 <- c()
length2 <- c()
for (i in row.names(non))
{print(i)
  length1<-rbind(length1,s[non[i,1],]$width)
  length2<-rbind(length2,s[non[i,2],]$width)
}
non$length1 <- length1
non$length2 <- length2
## adding overlap length
overlap<- c()
for (i in 1:nrow(non))
{ print(i)
  a1<- s[non[i,2],]$start 
  a2 <-  s[non[i,2],]$end 
  b1 <-  s[non[i,1],]$start
  b2 <- s[non[i,1],]$end 
  if (a2<=b2 & a1>=b1 & a2>b1) {overlap <- rbind(overlap, data.frame(id1=non[i,1],id2=non[i,2],overlapc="A",overlap=a2-a1))}
  if (a2>=b2 & a1<=b1 & b2>a1) {overlap<- rbind(overlap, data.frame(id1=non[i,1],id2=non[i,2],overlapc="B",overlap=b2-b1))}
  if (a2>b2 & a1<b2 & b1<a1) {overlap <- rbind(overlap, data.frame(id1=non[i,1],id2=non[i,2],overlapc="C",overlap=b2-a1))}
  if (a2<b2 & a1<b1 & b1<a2) {overlap <- rbind(overlap, data.frame(id1=non[i,1],id2=non[i,2],overlapc="D",overlap=a2-b1))}
}
dim(overlap)
dim(unique(overlap))
## list the repetitive gene lists and define 95% coverage as class E
non[which(non$length1==non$length2),]
#  V1        V2 promoterdis overlap length1 length2
#440  AT1G05560 AT1G05562        2135      NA    2136    2136
#989  AT3G21780 AT3G21781        1615      NA    1616    1616
# 1087 AT3G24330 AT3G24332        1502      NA    1503    1503
## keep three lines unique
overlap[overlap$id1=="AT1G05560",]
overlap <-overlap[-1154,]
overlap[overlap$id1=="AT1G05560",]$overlapc <- "E"
overlap[overlap$id1=="AT3G21780",]
overlap <-overlap[-975,]
overlap[overlap$id1=="AT3G21780",]$overlapc <- "E"
overlap[overlap$id1=="AT3G24330",]
overlap <-overlap[-927,]
overlap[overlap$id1=="AT3G24330",]$overlapc <- "E"
## successfully removed extra lines
dim(overlap)
non$overlap <- overlap$overlap
non$goverlap <- overlap$overlapc
## extract expression value from single RNA datas
d=row.names(sin1)
non$sin1<- ifelse(non$id2 %in%d & non$id1 %in%d, "1","0")
d=row.names(sin2)
non$sin2<- ifelse(non$id2 %in%d & non$id1 %in%d, "1","0")
d=row.names(sin3)
non$sin3<- ifelse(non$id2 %in%d & non$id1 %in%d, "1","0")
d=row.names(snuc)
non$snuc<- ifelse(non$id2 %in%d & non$id1 %in%d, "1","0")
non$sin1 <- as.integer(non$sin1)
non$sin2 <- as.integer(non$sin2)
non$sin3 <- as.integer(non$sin3)
non$snuc <- as.integer(non$snuc)
non <- non[!rowSums(non[,8:11])=="0",]
write.table(non, file="candidates_RNAPII_collision.txt", sep="\t", quote=F, row.names=F, col.names=T)
z <- non[non$snuc==1,]
m <- non[rowSums(non[,8:11])>=3,]
m <- m[order(m$overlap,decreasing = T),]
m <- m[order(m$promoterdis,decreasing = T),]
 s1 <- non[!non$snuc==0,]
 k1 <- c()
 k2 <- c()
 for (i in 1:nrow(s1)){
   print(i)
if(s1[i,1]%in%d & s1[i,2]%in%d) 
{k1 <-rbind(k1,snuc[s1[i,1],])
k2 <-rbind(k2,snuc[s1[i,2],])
}
 }
k1 <- k1[,-1]
k1$median <- apply(k1, 1, median)
head(k1)
k2<- k2[,-1]
k2$median <- apply(k2, 1, median)
s1$sex <- k1$median
s1$asex <- k2$median
##make sure median value is correct
k1["AT2G26360",]
s1 <- s1[order(s1$asex,decreasing = T),]
k1["AT4G31250",]
## AT5G66560, AT5G66558
## also put bulk expression 
dim(s1[s1$sex>0.1&s1$asex>0.1,])
dim(s1[s1$sex>0.05&s1$asex>0.05,])
write.table(s1, file="candidates_RNAPII_collision_snuc.txt", sep="\t", quote=F, row.names=F, col.names=T)

## for coding_anticoding
co_co <- read.table("/Users/xung0001/Desktop/coding_anticoding.txt",head=T,sep = ",")
 non <- co_co
colnames(non)[1:2] <- c("id1","id2")
row.names(mRNA) <- mRNA$gene_id
length <- c()
for (i in 1:nrow(non))
{print(i)
  if(mRNA[non[i,2],]$strand=="-")
  {length <-rbind(length,data.frame(id1=non[i,1],id2=non[i,2], length=mRNA[non[i,2],]$end-mRNA[non[i,1],]$start))}
  if(mRNA[non[i,2],]$strand=="+") {length <-rbind(length,data.frame(id1=non[i,1],id2=non[i,2], length=mRNA[non[i,1],]$end-mRNA[non[i,2],]$start))
  }
  
}
non$promoterdis <- length$length
non <- non[order(non$promoterdis,decreasing = F),]


length1 <- c()
length2 <- c()
for (i in row.names(non))
{print(i)
  length1<-rbind(length1,mRNA[non[i,1],]$width)
  length2<-rbind(length2,mRNA[non[i,2],]$width)
}
non$length1 <- length1
non$length2 <- length2  
overlap<- c()
for (i in 1:nrow(non))
{ print(i)
  a1 <-  mRNA[non[i,2],]$end 
  b1 <-  mRNA[non[i,1],]$end 
  a2 <- mRNA[non[i,2],]$start 
  b2 <- mRNA[non[i,1],]$start
  if (a1>=b1 & a2>=b2) {overlap <- rbind(overlap, data.frame(id1=non[i,1],id2=non[i,2],overlapc="A",overlap=b1-a2))}
  if (b1>=a1 & b2>=a2 ) {overlap<- rbind(overlap, data.frame(id1=non[i,1],id2=non[i,2],overlapc="A",overlap=a1-b2))}
  if (b1>=a1 &  b2<=a2) {overlap <- rbind(overlap, data.frame(id1=non[i,1],id2=non[i,2],overlapc="B",overlap=a1-a2))}
  if (a1>=b1 &  b2>=a2) {overlap <- rbind(overlap, data.frame(id1=non[i,1],id2=non[i,2],overlapc="B",overlap=b1-b2))}
}
dim(overlap)
overlap <- unique(overlap)
row.names(overlap) <- overlap$id1 
overlap <- overlap[-535,]
overlap <- overlap[-1098,]
overlap <- overlap[-389,]
overlap <- overlap[-436,]
overlap <- overlap[-728,]
overlap <- overlap[-1402,]
non$overlap <- overlap$overlap
non$overlapc <- overlap$overlapc
## extract expression value from single RNA datas
d=row.names(sin1)
non$sin1<- ifelse(non$id2 %in%d & non$id1 %in%d, "1","0")
d=row.names(sin2)
non$sin2<- ifelse(non$id2 %in%d & non$id1 %in%d, "1","0")
d=row.names(sin3)
non$sin3<- ifelse(non$id2 %in%d & non$id1 %in%d, "1","0")
d=row.names(snuc)
non$snuc<- ifelse(non$id2 %in%d & non$id1 %in%d, "1","0")
non$sin1 <- as.integer(non$sin1)
non$sin2 <- as.integer(non$sin2)
non$sin3 <- as.integer(non$sin3)
non$snuc <- as.integer(non$snuc)
non <- non[!rowSums(non[,8:11])=="0",]
write.table(non, file="candidates_RNAPII_collision_c_coding.txt", sep="\t", quote=F, row.names=F, col.names=T)
z <- non[non$snuc==1,]
m <- non[rowSums(non[,8:11])>=3,]
m <- m[order(m$overlap,decreasing = T),]
m <- m[order(m$promoterdis,decreasing = T),]
# 923 pairs in single nucleus 
s1 <- non[!non$snuc==0,]
k1 <- c()
k2 <- c()
for (i in 1:nrow(s1)){
  print(i)
  if(s1[i,1]%in%d & s1[i,2]%in%d) 
  {k1 <-rbind(k1,snuc[s1[i,1],])
  k2 <-rbind(k2,snuc[s1[i,2],])
  }
}
k1 <- k1[,-1]
k1$median <- apply(k1[,-1], 1, median)
head(k1)
k2<- k2[,-1]
k2$median <- apply(k2[,-1], 1, median)
s1$sex <- k1$median
s1$asex <- k2$median
##make sure median value is correct
k1["AT3G58520",]
k2["AT3G58530",]
## AT5G66560, AT5G66558
# 923 pairs 
s1 <- s1[order(s1$asex,decreasing = T),]
## also put bulk expression 
write.table(s1, file="candidates_RNAPII_collision_snuc_coding.txt", sep="\t", quote=F, row.names=F, col.names=T)
s1 <- s1[order(s1$overlap,decreasing = T),]
s2 <- s1[s1$overlap>500,] #188 paris 
s3 <- read.table("candidates_RNAPII_collision_snuc.txt",head=T,sep = "\t")
# coding_noncoding has more pairs with bigger overlap distance 
s4 <- s3[s3$overlap>500,] #234 pairs
s2<- s2[order(s2$asex,decreasing = T),]
s4<- s4[order(s4$overlap,decreasing = T),]
dim(s1[rowSums(s1[,12:13])>0,])