
#1.read files
setwd("C:\\Users\\lenovo\\Desktop\\survival analysis\\cluster")
Cluster=read.table("cluster.txt",header = T,row.names = 1,check.names = F)
colnames(Cluster)=c("cluster")

time=read.table("time.txt",header = T,row.names = 1,check.names = F)
time$futime=time$futime/365

#2.insect two files 
same_id=intersect(row.names(Cluster),row.names(time))
time=time[same_id,]
clu_time=cbind(time,Cluster)

#3.survival analysis 
library(survival)
diff=survdiff(Surv(as.numeric(futime,fustat))~cluster,data=clu_time)
length=length(levels(factor(clu_time$cluster)))

pValue=1-pchisq(diff$chisq, df=length-1)

if(pValue<0.001){
  pValue="p<0.001"
}else{
  pValue=paste0("p=",sprintf("%.03f",pValue))
}

fit <- survfit(Surv(as.numeric(futime, fustat)) ~ cluster, data = clu_time)

#4.paint the survival curve 

library(dplyr)
library(survminer)

pdf(file="survival.pdf",onefile = FALSE,width=5.5,height=6.5)
ggsurvplot(fit,pval = T,
           conf.int = T,
           conf.int.style="ribbon",
           surv.median.line = "hv",
           conf.int.alpha=0.1,
           palette = "jco",
           risk.table = T)
dev.off()