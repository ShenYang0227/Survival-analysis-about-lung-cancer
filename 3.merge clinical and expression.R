##merge the clinical information with the expression of m6A_lncRNA
#1.read files
setwd("C:\\Users\\lenovo\\Desktop\\gdc_download_20230510_123138.115629")
m6aLnc=read.table("m6aLncRNA.txt",header = TRUE,check.names = FALSE)
row.names(m6aLnc)=make.unique(m6aLnc$id)
m6aLnc=m6aLnc[,-1]

library(limma)
m6aLnc=avereps(m6aLnc)

clin=read.table("time.txt",header = TRUE,row.names=1,check.names = FALSE)
clin$futime=clin$futime/365

#2.find out the same id between m6aLnc and clin
colnames(m6aLnc)=substring(colnames(m6aLnc),1,12)

same_id=intersect(colnames(m6aLnc),row.names(clin))
m6aLnc_exp=m6aLnc[,same_id]
m6aLnc_exp=t(m6aLnc_exp)
clin_exp=clin[same_id,]
m6aLnc_clin=cbind(id=row.names(m6aLnc_exp),clin_exp,m6aLnc_exp)

#3.save and write
write.table(m6aLnc_clin,file = "m6alnc_time.txt",sep="\t",quote = F,row.names = F)