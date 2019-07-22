source("https://raw.githubusercontent.com/Shicheng-Guo/GscRbasement/master/GscTools.R")
library("metafor")

setwd("~/hpc/methylation/Pancancer/RNA-seq")
load("rnaseqdata.pancancer.RData")

TCGAProjects=c("BLCA","BRCA","CESC","CHOL","COAD","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LIHC","LUAD","LUSC","PAAD","PCPG","PRAD","READ","SARC","STAD","THCA","THYM","UCEC")
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

idx<-which(phen$phen2==1 | phen$phen2==11)
phen<-phen[idx,]
input<-rnaseqdata[,idx]

idx<-which(phen$pid %in% paste("TCGA-",TCGAProjects,sep=""))
phen<-phen[idx,]
input<-input[,idx]
input[1:5,1:5]

input<-log(input+1,2)
input<-RawNARemove(input)
input<-RawZeroRemove(input)

Seq<-paste(phen$project_id,phen$phen2,sep="-")
rlt<-c()
coll<-c()
i=match("ENSG00000213754.2",rownames(input))
i=match("ENSG00000231246.1",rownames(input))
panc<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/PANC/master/extdata/panc.txt",head=T)
ii<-unlist(lapply(rownames(panc),function(x) grep(x,rownames(input))))

for(i in ii){
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
  Source<-unlist(lapply(strsplit(names(m1i),"-"),function(x) x[2]))
  output<-data.frame(cbind(n1i,m1i,sd1i,n2i,m2i,sd2i))
  output$source=Source
  output<-na.omit(output)
  es<-escalc(m1i=m1i, sd1i=sd1i, n1i=n1i, m2i=m2i, sd2i=sd2i, n2i=n2i,measure="MD",data=output)
  md <- rma(es,slab=source,method = "REML", measure = "SMD",data=output)
  rlt<-rbind(rlt,c(i,md$beta,md$pval,md$ci.lb,md$ci.ub,md$I2,md$tau2))
  coll<-c(coll,i)
  m<-metagen(yi,seTE=vi,data = es,
             comb.fixed = TRUE,
             comb.random = TRUE,
             prediction=F,
             sm="SMD")
  pdf(paste("./OS_HR/",rownames(input)[i],".SMD.PANC.pdf",sep=""))
  forest(m,leftlabs = rownames(HR),
         lab.e = "Intervention",
         pooled.totals = FALSE,
         smlab = "",studlab=rownames(HR),
         text.random = "Overall effect",
         print.tau2 = FALSE,
         col.diamond = "blue",
         col.diamond.lines = "black",
         col.predict = "red",
         print.I2.ci = TRUE,
         digits.sd = 2,fontsize=8)
  dev.off()
}

library("meta")
library("metafor")
library("survival")
library("survminer")
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
panc<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/PANC/master/extdata/panc.txt",head=T)
ii<-unlist(lapply(rownames(panc),function(x) grep(x,rownames(input))))

i=match("ENSG00000213754.2",rownames(input))
i=match("ENSG00000231246.1",rownames(input))

for(i in ii){
  Z<-c()
  for(z in quantile(input[i,],seq(0, 1, 0.025))[2:40]){
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
    print(c(z,P))
    m<-metagen(HR[,1],seTE=HR[,3],comb.fixed = TRUE,comb.random = TRUE,prediction=F,sm="HR")
    Z<-rbind(Z,c(z,m$pval.fixed))
    print(c(z,m$pval.fixed))
  }
  
  thres<-Z[which.min(Z[,2]),1]
  HR<-c()
  for(TCGAProject in TCGAProjects){
    newdata<-input[,phen$project_id==paste("TCGA-",TCGAProject,sep="")]
    xphen<-phen[phen$project_id==paste("TCGA-",TCGAProject,sep=""),]
    dat<-data.frame(Rna=newdata[i,],xphen)
    dat$Rna[dat$Rna<=thres]<-0
    dat$Rna[dat$Rna>thres]<-1
    hr.fit<-summary(coxph(Surv(week,censored)~Rna,dat))
    hr1=hr.fit$coefficients[1,]
    hr2=hr.fit$conf.int[1,]
    HR<-rbind(HR,c(hr1,hr2[3],hr2[4]))
    
    fit <- survfit(Surv(week,censored)~Rna, data = dat)
    survp<-ggsurvplot(fit, data = dat,conf.int = F,pval = TRUE,
                      fun = "pct",risk.table = TRUE,size = 1,linetype = "strata",
                      palette = c("#E7B800","#2E9FDF"),
                      legend = "bottom",legend.title = rownames(input)[i],
                      legend.labs = c("High-expression","Low-expression"))
    ggsave(file = paste("./OS_HR/",rownames(input)[i],"_",TCGAProject,"_KM.pdf",sep=""), survp$plot)
    
  }
  print(i)
  rownames(HR)<-TCGAProjects

  m<-metagen(HR[,1],seTE=HR[,3],
             comb.fixed = TRUE,
             comb.random = TRUE,
             prediction=F,
             sm="HR")
  
  write.table(m,file=paste("./OS_HR/",rownames(input)[i],".OS.HR.m.PANC.txt",sep=""),sep="\t",col.names = NA,row.names = T,quote=F)
  write.table(HR,file=paste("./OS_HR/",rownames(input)[i],".OS.HR.n.PANC.txt",sep=""),sep="\t",col.names = NA,row.names = T,quote=F)
  pdf(paste("./OS_HR/",rownames(input)[i],".OS.HR.PANC.pdf",sep=""))
  forest(m,leftlabs = rownames(HR),
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
  
}


file<-list.files(pattern="*.txt")
for(i in file){
  data<-na.omit(read.table(i,sep="\t",head=T))
  if(nrow(data)>1){
  data <- data[!is.infinite(rowSums(data[,6:7])),]
  fit<-t.test((1-data[,6])*(1-data[,7]))
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
write.table(HR,file=paste("./OS_HR/",rownames(input)[i],".OS.HR.Pancancer.txt",sep=""),sep="\t",col.names = NA,row.names = T,quote=F)

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


i=match("ENSG00000231246.1",rownames(input))

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
  m<-metagen(HR[,1],seTE=HR[,3],comb.fixed = TRUE,comb.random = TRUE,prediction=F,sm="HR")
  print(c(z,m$pval.fixed))
  
}

HR<-c()
z=5000
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
write.table(HR,file=paste("./OS_HR/",rownames(input)[i],".OS.HR.Pancancer.txt",sep=""),sep="\t",col.names = NA,row.names = T,quote=F)

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

pdf("./OS_HR/ENSG00000231246.1.pdf")
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










