setwd("/mnt/bigdata/Genetic/Projects/shg047/methylation/Pancancer/RNA-seq")

Symbol2ENSG<-function(Symbol){
  db<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/AnnotationDatabase/master/ENSG.ENST.ENSP.Symbol.hg19.bed",sep="\t")
  ENSG<-as.character(db[match(Symbol,db$V4),8])
  ENSG<-na.omit(data.frame(Symbol,ENSG))
  return(ENSG)
}
ENSG2Symbol<-function(ENSG){
  db<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/AnnotationDatabase/master/ENSG.ENST.ENSP.Symbol.hg19.bed",sep="\t")
  ENSG<-unlist(lapply(strsplit(ENSG,split="[.]"),function(x) x[1]))
  Symbol<-db[match(as.character(ENSG),db$V8),4]
  return(Symbol)
}


ENSG<-Symbol2ENSG(c("TP53","KRAS","TTN","MUC16","CSMD3","PIK3CA","CDK2"))
xgene<-c("ENSG00000213754",as.character(ENSG[,2]))
xsel<-unlist(lapply(xgene,function(x) grep(x,rownames(input))))
name<-as.character(ENSG2Symbol(xgene))
name<-c("PANC754","TP53","KRAS","TTN","MUC16","CSMD3","PIK3CA","CDK2")
temp<-input[xsel,]
rownames(temp)<-name
value<-cor(t(temp))
value
pdf("../../PANC754.TP53.heatmap.pdf")
HeatMap(cor(t(temp)),cexRow=0.75,cexCol=0.75)
dev.off()
write.table(value,file="../../PANC754.TP53.correlation.txt",sep="\t",quote=F,col.names = NA,row.names = T)
P<-c()
corvalue<-c()
for(i in 2:nrow(temp)){
  xx<-cor.test(temp[1,],temp[i,])
  corvalue<-c(corvalue,xx$estimate)
  P<-c(P,xx$p.value)
}

rlt<-data.frame(corvalue,P)
rownames(rlt)<-rownames(temp)[2:nrow(temp)]
write.table(rlt,file="../../PANC754.TP53.co-expression.pvalue.txt",sep="\t",col.names = NA,row.names = T,quote=F)
