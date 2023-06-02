#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("ggplot2")
#install.packages("ggpubr")


library(limma)
library(ggplot2)
library(ggpubr)
expFile="symbol.txt"            
clusterFile="cluster.txt"       
gene="CD274"                    
showName="PD-L1"                
setwd("C:\\Users\\lenovo\\Desktop\\15 geneDiff")     


rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)
data=t(data[gene,,drop=F])
data=log2(data+1)


group=sapply(strsplit(rownames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
conNum=length(group[group==1])       
treatNum=length(group[group==0])    
Type=c(rep(1,conNum), rep(2,treatNum))



exp=cbind(data, Type)
exp=as.data.frame(exp)
colnames(exp)=c("gene", "Type")
exp$Type=ifelse(exp$Type==1, "Normal", "Tumor")
#exp$gene[exp$gene>quantile(exp$gene,0.99)]=quantile(exp$gene,0.99)


group=levels(factor(exp$Type))
exp$Type=factor(exp$Type, levels=group)
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}


bioCol=c("#0066FF","#FF0000","#FF9900","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(levels(factor(exp$Type)))]

boxplot=ggboxplot(exp, x="Type", y="gene", color="Type",
		          xlab="",
		          ylab=paste0(showName, " log2(expression)"),
		          legend.title="Type",
		          palette = bioCol,
		          add = "jitter")+ 
	stat_compare_means(comparisons=my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")

pdf(file=paste0(gene,".diff.pdf"), width=5, height=4.5)
print(boxplot)
dev.off()



cluster=read.table(clusterFile, header=F, sep="\t", check.names=F, row.names=1)
sameSample=intersect(row.names(data), row.names(cluster))
data=data.frame(gene=data[sameSample,], Cluster=cluster[sameSample,])
data$Cluster=paste0("Cluster", data$Cluster)
#data$gene[data$gene>quantile(data$gene,0.99)]=quantile(data$gene,0.99)


group=levels(factor(data$Cluster))
data$Cluster=factor(data$Cluster, levels=group)
comp=combn(group, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}


bioCol=c("#0066FF","#FF0000","#FF9900","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(levels(factor(data$Cluster)))]

boxplot=ggboxplot(data, x="Cluster", y="gene", color="Cluster",
		          xlab="",
		          ylab=paste0(showName, " log2(expression)"),
		          legend.title="Cluster",
		          palette = bioCol,
		          add = "jitter")+ 
	stat_compare_means(comparisons=my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")

pdf(file=paste0(gene,".Cluster.pdf"), width=5, height=4.5)
print(boxplot)
dev.off()


