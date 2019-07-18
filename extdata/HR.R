install.packages("metafor")
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

for(i in ii){
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
  }
  print(i)
  rownames(HR)<-TCGAProjects
  write.table(HR,file=paste("./OS_HR/",rownames(input)[i],".OS.HR.PANC.txt",sep=""),sep="\t",col.names = NA,row.names = T,quote=F)
  print(HR)
  
  m<-metagen(HR[,1],seTE=HR[,3],
             comb.fixed = TRUE,
             comb.random = TRUE,
             prediction=F,
             sm="HR")
  
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
