#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("corrplot")


#ÒýÓÃ°ü
library(limma)
library(corrplot)
expFile="symbol.txt"        
lncFile="uniSigExp.txt"     
gene="CD274"                
showName="PD-L1"            
setwd("C:\\Users\\lenovo\\Desktop\\16 geneCor")     
 
#1.read files
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)



lncRT=read.table(lncFile, header=T, sep="\t", check.names=F, row.names=1)
lncRNA=c(gene, colnames(lncRT)[3:ncol(lncRT)])
sameGene=intersect(lncRNA, row.names(data))
data=data[sameGene,]
row.names(data)[1]=showName

#2. make correction matrix 
data=t(data)
M=cor(data)
res1=cor.mtest(data, conf.level = 0.95)

#3. paint correction map 
pdf(file="cor.pdf", width=8, height=8)
corrplot(M,
         order="original",
         method = "circle",
         type = "upper",
         tl.cex=0.8, pch=T,
         p.mat = res1$p,
         insig = "label_sig",
         pch.cex = 1.6,
         sig.level=0.05,
         number.cex = 1,
         col=colorRampPalette(c("blue", "white", "red"))(50),
         tl.col="black")
dev.off()


