#install.packages("pheatmap")


library(pheatmap)               
ClusterFile="cluster.txt"       
cliFile="clinical.txt"          
expFile="uniSigExp.txt"         
setwd("C:\\Users\\lenovo\\Desktop\\14 clusterCli")

#1. read related files 
Cluster=read.table(ClusterFile, header=F, sep="\t", check.names=F, row.names=1)
colnames(Cluster)=c("Cluster")
Cluster=Cluster[order(Cluster$Cluster),,drop=F] 


cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)

#2.intersect two files 
samSample=intersect(row.names(Cluster), row.names(cli))
Cluster=Cluster[samSample,"Cluster",drop=F]
cli=cli[samSample,,drop=F]
Type=cbind(Cluster, cli)
Type$Cluster=paste0("Cluster", Type$Cluster)

#3. correction analysis about clinical informations 
sigVec=c("Cluster")
for(clinical in colnames(Type[,2:ncol(Type)])){
	data=Type[c("Cluster", clinical)]
	colnames(data)=c("Cluster", "clinical")
	data=data[(data[,"clinical"]!="unknow"),]
	tableStat=table(data)
	stat=chisq.test(tableStat)
	pvalue=stat$p.value
	Sig=ifelse(pvalue<0.001,"***",ifelse(pvalue<0.01,"**",ifelse(pvalue<0.05,"*","")))
	sigVec=c(sigVec, paste0(clinical, Sig))
	#print(paste(clinical, pvalue, Sig, sep="\t"))
}
colnames(Type)=sigVec

#4. obtain the expression 
exp=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
data=t(exp[,(3:ncol(exp))])
data=data[,row.names(Type)]

#5. heatmap 
colorList=list()
#Type=Type[apply(Type,1,function(x)any(is.na(match('unknow',x)))),,drop=F]
bioCol=c("#0066FF","#FF0000","#FF9900","#ed1299", "#0dbc21", "#246b93", "#cc8e12", "#d561dd", "#c93f00", 
         "#ce2523", "#f7aa5d", "#9ed84e", "#39ba30", "#6ad157", "#373bbf", "#a1ce4c", "#ef3bb6", "#d66551",
         "#1a918f", "#7149af", "#ff66fc", "#2927c4", "#57e559" ,"#8e3af4" ,"#f9a270" ,"#22547f", "#db5e92",
         "#4aef7b", "#e86502",  "#99db27", "#e07233", "#8249aa","#cebb10", "#03827f", "#931635", "#ff523f",
         "#edd05e", "#6f25e8", "#0dbc21", "#167275", "#280f7a", "#6373ed", "#5b910f" ,"#7b34c1" ,"#0cf29a" ,"#d80fc1",
         "#dd27ce", "#07a301", "#ddd53e",  "#391c82", "#2baeb5","#925bea", "#09f9f5",  "#63ff4f")
j=0
for(cli in colnames(Type[,1:ncol(Type)])){
	cliLength=length(levels(factor(Type[,cli])))
	cliCol=bioCol[(j+1):(j+cliLength)]
	j=j+cliLength
	names(cliCol)=levels(factor(Type[,cli]))
	if("unknow" %in% levels(factor(Type[,cli]))){
		cliCol["unknow"]="grey75"}
	colorList[[cli]]=cliCol
}


data=log2(data+1)
pdf("heatmap.pdf", width=7.5, height=5)
pheatmap(data,
         annotation=Type,
         annotation_colors = colorList,
         color = colorRampPalette(c(rep("blue",5), "white", rep("red",5)))(100),
         cluster_cols =F,
         cluster_rows =T,
         scale="row",
         show_colnames=F,
         show_rownames=T,
         fontsize=6,
         fontsize_row=6,
         fontsize_col=6)
dev.off()


