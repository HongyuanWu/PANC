
cd ~/hpc/tcga/PancRNA
data<-read.table("TCGA-SMD-DGE-Meta-Pvalue-2019.txt",head=T,sep="\t")
new<-data[which(is.na(data[,ncol(data)])),]

new[grep("ENSG00000213754",new[,1]),]

ENSG2Symbol<-function(ENSG){
  db<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/AnnotationDatabase/master/ENSG.ENST.ENSP.Symbol.hg19.bed",sep="\t",as.is=T)
  ENSG<-unlist(lapply(strsplit(ENSG,split="[.]"),function(x) x[1]))
  symbol<-db[db$V8 %in% as.character(ENSG),]
  return(symbol)
}


ensg2bed<-function(ENSG){
  db<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/AnnotationDatabase/master/hg19/ENSG.ENST.hg19.txt",as.is=T,head=F)
  ENSG<-unlist(lapply(strsplit(ENSG,split="[.]"),function(x) x[1]))
  bed<-unique(db[db$V5 %in% as.character(ENSG),c(1,2,3,5)])
  return(bed)
}


bed<-ensg2bed(as.character(new[,1]))
xsel<-unlist(apply(bed,1,function(x) grep(x[4],new[,1])))
output<-data.frame(bed,new[xsel,])
library("CMplot")
cminput<-data.frame(SNP=output$V5,Chromosome=output$V1,Position=output$V2,trait1=output[,13])
CMplot(cminput,plot.type="b",ylim=20,LOG10=TRUE,threshold=NULL,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE,width=14,height=6)

cminput[grep("ENSG00000213754",cminput[,1]),]
write.table(new,file="pancancer.meta.dge.11529.txt",sep="\t",quote=F,row.name=T,col.names=NA)
