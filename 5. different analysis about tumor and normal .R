#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("pheatmap")
#install.packages("reshape2")
#install.packages("ggpubr")


#引用包
library(limma)
library(pheatmap)
library(reshape2)
library(ggpubr)
lncFile="uniSigExp.txt"     #预后相关的lncRNA列表
expFile="lncRNA.txt"        #lncRNA表达文件
setwd("C:\\Users\\lenovo\\Desktop\\11 预后相关lncRNA的热图")       #设置工作目录

#读取输入文件
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

#提取预后lncRNA表达量
lncRNA=read.table(lncFile, header=T, sep="\t", check.names=F, row.names=1)
data=data[colnames(lncRNA)[3:ncol(lncRNA)],]
data=log2(data+1)
exp=data

#正常和肿瘤数目
group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
conNum=length(group[group==1])       #正常组样品数目
treatNum=length(group[group==0])     #肿瘤组样品数目
sampleType=c(rep(1,conNum), rep(2,treatNum))

#基因差异分析
sigVec=c()
for(i in row.names(data)){
	test=t.test(data[i,] ~ sampleType)
	pvalue=test$p.value
	Sig=ifelse(pvalue<0.001,"***",ifelse(pvalue<0.01,"**",ifelse(pvalue<0.05,"*","")))
	sigVec=c(sigVec, paste0(i, Sig))
}
row.names(data)=sigVec


#热图可视化
Type=c(rep("Normal",conNum), rep("Tumor",treatNum))
names(Type)=colnames(data)
Type=as.data.frame(Type)
data=log2(data+1)
pdf("heatmap.pdf", width=7.5, height=4.7)
pheatmap(data,
         annotation=Type,
         color = colorRampPalette(c(rep("blue",5), "white", rep("red",5)))(100),
         cluster_cols =F,
         cluster_rows =T,
         scale="row",
         show_colnames=F,
         show_rownames=T,
         fontsize=6,
         fontsize_row=7,
         fontsize_col=6)
dev.off()

#把数据转换成ggplot2输入文件
exp=as.data.frame(t(exp))
exp=cbind(exp, Type=sampleType)
exp$Type=ifelse(exp$Type==1, "Normal", "Tumor")
data=melt(exp, id.vars=c("Type"))
colnames(data)=c("Type", "Gene", "Expression")

#绘制箱线图
p=ggboxplot(data, x="Gene", y="Expression", color = "Type", 
	     ylab="log2(Gene expression)",
	     xlab="",
	     legend.title="Type",
	     palette = c("blue", "red"),
	     width=1)
p=p+rotate_x_text(60)
p1=p+stat_compare_means(aes(group=Type),
	      method="wilcox.test",
	      symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
	      label = "p.signif")

#输出箱线图
pdf(file="boxplot.pdf", width=15, height=8)
print(p1)
dev.off()


######生信自学网: https://www.biowolf.cn/
######课程链接1: https://shop119322454.taobao.com
######课程链接2: https://ke.biowolf.cn
######课程链接3: https://ke.biowolf.cn/mobile
######光俊老师邮箱：seqbio@foxmail.com
######光俊老师微信: eduBio
