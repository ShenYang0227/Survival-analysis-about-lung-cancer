##get the m6A-LncRNA by co-expression analysis
#1.read related files
library(limma)
setwd("C:\\Users\\lenovo\\Desktop\\gdc_download_20230510_123138.115629")
lncRNA=read.table("lncRNA.txt",header = TRUE,check.names = F,row.names = 1)
lncRNA=avereps(lncRNA)

m6A=read.table("m6A_matrix.txt",header = TRUE,check.names = F,row.names = 1)
m6A=avereps(m6A)

#2.grouping
table(substring(colnames(lncRNA),14,15)=="11")
group=c(rep("tumor",151),rep("normal",9))

#3.correction analysis
corFilter=0.25
pFilter=0.05
out=data.frame()
for (i in row.names(lncRNA)){
  if (sd(lncRNA[i,])>0.1){
    test=wilcox.test(as.numeric(lncRNA[i,])~group)
    if (test$p.value<0.05){
      for (j in row.names(m6A)){
        x=as.numeric(lncRNA[i,])
        y=as.numeric(m6A[j,])
        corT=cor.test(x,y)
        cor=corT$estimate
        pvalue=corT$p.value
        if ((cor>corFilter)&(pvalue<pFilter)){
          out=rbind(out,cbind(m6A=j,lncRNA=i,cor,pvalue,regulation="positive"))
        }
        if ((cor<corFilter)&(pvalue<pFilter)){
          out=rbind(out,cbind(m6A=j,lncRNA=i,cor,pvalue,regulation="negative"))
        }
      }
    }
    
  }
}

#4.find the significant genes and extract the expression about them 
m6aLncRNA=out[,"lncRNA"]
m6aLncRNAexp=lncRNA[m6aLncRNA,]
m6aLncRNAexp=cbind(id=row.names(m6aLncRNAexp),m6aLncRNAexp)

#5. save results 
write.table(out,file = "corTest.txt",sep = "\t",row.names = F,quote = F)
write.table(m6aLncRNAexp,file = "m6aLncRNA.txt",row.names = F,quote = F)