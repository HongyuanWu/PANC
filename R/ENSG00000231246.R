
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



