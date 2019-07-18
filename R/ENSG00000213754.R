source("https://raw.githubusercontent.com/Shicheng-Guo/GscRbasement/master/GscTools.R")
load("rnaseqdata.pancancer.RData")
input<-log(input+1,2)
input<-RawNARemove(input)
input<-RawZeroRemove(input)
Seq<-paste(phen$project_id,phen$bin,sep="-")
rlt<-c()
coll<-c()
for(i in 1:nrow(input)){
  print(i)
  mean<-tapply(as.numeric(input[i,]),Seq,function(x) mean(x,na.rm=T))
  sd<-tapply(as.numeric(input[i,]),Seq,function(x) sd(x,na.rm=T))
  num<-tapply(as.numeric(input[i,]),Seq,function(x) length(x))
  m1i=mean[seq(1,length(mean),by=2)]
  m2i=mean[seq(2,length(mean),by=2)]
  sd1i=sd[seq(1,length(mean),by=2)]
  sd2i=sd[seq(2,length(mean),by=2)]
  n1i=num[seq(1,length(mean),by=2)]
  n2i=num[seq(2,length(mean),by=2)]
  Source<-unlist(lapply(strsplit(names(m1i),"_"),function(x) x[1]))
  output<-data.frame(cbind(n1i,m1i,sd1i,n2i,m2i,sd2i))
  output$source=Source
  output<-na.omit(output)
  es<-escalc(m1i=m1i, sd1i=sd1i, n1i=n1i, m2i=m2i, sd2i=sd2i, n2i=n2i,measure="MD",data=output)
  md <- rma(es,slab=source,method = "REML", measure = "SMD",data=output)
  rlt<-rbind(rlt,c(i,md$beta,md$pval,md$ci.lb,md$ci.ub,md$I2,md$tau2))
  coll<-c(coll,i)
}
rownames(rlt)<-rownames(input)[coll]
colnames(rlt)<-c("idx","beta","pval","cilb","ciub","i2","tau2")
rlt<-data.frame(rlt)



library("metafor")
library("survival")
library("survminer")
install.packages("metafor")
setwd("~/hpc/methylation/Pancancer/RNA-seq")
load("rnaseqdata.pancancer.RData")

TCGAProjects=c("BLCA","BRCA","CESC","CHOL","COAD","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LIHC","LUAD","LUSC","PAAD","PCPG","PRAD","READ","SARC","STAD","THCA","THYM","UCEC")
OS<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/TCGA/OverallSurvivalTime.txt",head=T,sep="\t")
phen1=read.table("~/hpc/methylation/TCGA-clinical-11093.tsv",header = T,sep="\t")
phen2=read.table("~/hpc/methylation/Pancancer/RNA-seq/File_metadata2.txt",header = T,sep="\t")
head(phen1)
head(phen2)
colnames(rnaseqdata)<-unlist(lapply(strsplit(colnames(rnaseqdata),"/"),function(x) x[2]))
phen<-data.frame(phen2,phen1[match(phen2$cases.0.case_id,phen1$case_id),])
phen$file_name=gsub(".gz","",phen$file_name)

# prepare phenotype information
phen<-phen[match(colnames(rnaseqdata),phen$file_name),]
phen$phen4<-id2phen4(phen$cases.0.samples.0.submitter_id)
phen$phen3<-id2phen3(phen$cases.0.samples.0.submitter_id)
phen$phen2<-id2bin(phen$cases.0.samples.0.submitter_id)
phen$pid<-phen$project_id
head(phen)
# enroll all the cancer samples and remove normal samples
idx<-which(c(phen$phen2==1))
phen<-phen[idx,]
input<-rnaseqdata[,idx]
input[1:5,1:5]

# match survival information
idx<-na.omit(match(OS$submitter_id,phen$phen3))
input<-input[,idx]
phen<-phen[idx,]
phen<-data.frame(phen,OS[match(phen$phen3,OS$submitter_id),])

phen$censored<-as.numeric(! phen$censored)
phen$week=phen$time/7

for(i in 1:nrow(input)){
  HR<-c()
  for(TCGAProject in TCGAProjects){
    newdata<-input[,phen$project_id==paste("TCGA-",TCGAProject,sep="")]
    xphen<-phen[phen$project_id==paste("TCGA-",TCGAProject,sep=""),]
    dat<-data.frame(Rna=newdata[i,],xphen)
    dat$Rna[dat$Rna<=median(dat$Rna)]<-0
    dat$Rna[dat$Rna>median(dat$Rna)]<-1
    hr.fit<-summary(coxph(Surv(week,censored)~Rna,dat))
    hr1=hr.fit$coefficients[1,]
    hr2=hr.fit$conf.int[1,]
    HR<-rbind(HR,c(hr1,hr2[3],hr2[4]))
  }
  print(i)
  rownames(HR)<-TCGAProjects
  write.table(HR,file=paste("./OS_HR/",rownames(input)[i],".OS.HR.Pancancer.txt",sep=""),sep="\t",col.names = NA,row.names = T,quote=F)
  print(HR)
}


file<-list.files(pattern="*.txt")
for(i in file){
  data<-na.omit(read.table(i,sep="\t",head=T))
  if(nrow(data)>1){
  data <- data[!is.infinite(rowSums(data[,7:8])),]
  fit<-t.test((1-data[,7])*(1-data[,8]))
  if(fit$p.value<10^-3){
  print(i)
    }
  }
}



BiocManager::install("survcomp")
BiocManager::install("forestplot")
require(rmeta)
require(survcomp)
require(forestplot)
library("meta")

i=match("ENSG00000213754.2",rownames(input))

Z<-c()
for(z in quantile(newdata[i,],seq(0, 1, 0.025))[2:40]){
  HR<-c()
  for(TCGAProject in TCGAProjects){
  newdata<-input[,phen$project_id==paste("TCGA-",TCGAProject,sep="")]
  xphen<-phen[phen$project_id==paste("TCGA-",TCGAProject,sep=""),]
  dat<-data.frame(Rna=newdata[i,],xphen)
  dat$Rna[dat$Rna<=z]<-0
  dat$Rna[dat$Rna>z]<-1
  hr.fit<-summary(coxph(Surv(week,censored)~Rna,dat))
  hr1=hr.fit$coefficients[1,]
  hr2=hr.fit$conf.int[1,]
  HR<-rbind(HR,c(hr1,hr2[3],hr2[4]))
}
  HR<-na.omit(HR)
  HR <- HR[!is.infinite(rowSums(HR[,6:7])),]
  P<-t.test((1-HR[,6])*(1-HR[,7]))$p.value
  Z<-rbind(Z,c(z,P))
  print(c(z,P))
}

HR<-c()
z=3200
for(TCGAProject in TCGAProjects){
  newdata<-input[,phen$project_id==paste("TCGA-",TCGAProject,sep="")]
  xphen<-phen[phen$project_id==paste("TCGA-",TCGAProject,sep=""),]
  dat<-data.frame(Rna=newdata[i,],xphen)
  dat$Rna[dat$Rna<=z]<-0
  dat$Rna[dat$Rna>z]<-1
  hr.fit<-summary(coxph(Surv(week,censored)~Rna,dat))
  hr1=hr.fit$coefficients[1,]
  hr2=hr.fit$conf.int[1,]
  HR<-rbind(HR,c(hr1,hr2[3],hr2[4]))
}
rownames(HR)<-TCGAProjects

library(forestplot)

clrs <- fpColors(box="royalblue",line="darkblue", summary="royalblue")
tabletext <- list(c(NA, rownames(HR)),append(list(expression(HR)), sprintf("%.3f", HR[,2])))
tabletext
pdf(paste("./OS_HR/",rownames(input)[i],".OS.HR.Pancancer.pdf",sep=""))
forestplot(tabletext, rbind(rep(NA, 3), HR[,c(2,6,7)]),col=clrs,xlab="EQ-5D index",zero =1)
dev.off()

m<-metagen(HR[,1],seTE=HR[,3],
           comb.fixed = TRUE,
           comb.random = TRUE,
           prediction=F,
           sm="HR")

pdf("./OS_HR/ENSG00000213754.2.forest.f11.pdf")
forest(m,leftlabs = rownames(HR),xlim=c(0.1,3),
       lab.e = "Intervention",
       pooled.totals = FALSE,
       smlab = "",studlab=rownames(HR),
       text.random = "Overall effect",
       print.tau2 = FALSE,
       col.diamond = "blue",
       col.diamond.lines = "black",
       col.predict = "red",
       print.I2.ci = TRUE,
       digits.sd = 2,fontsize=9)
dev.off()

