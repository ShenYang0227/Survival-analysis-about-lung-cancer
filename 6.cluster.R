
#1 read file 
setwd("C:\\Users\\lenovo\\Desktop\\survival analysis\\cluster")
exp=read.table("surSigExp.txt",header = T,row.names = 1,check.names = F)
exp=exp[,3:ncol(exp)]
exp=t(exp)


#2.cluster 
library(ConsensusClusterPlus)
maxK=9
result=ConsensusClusterPlus(exp,
                            maxK=maxK,
                            reps=50,
                            pItem=0.8,
                            pFeature=1,
                            title="C:\\Users\\lenovo\\Desktop\\survival analysis\\cluster",
                            clusterAlg="km",
                            distance="euclidean",
                            seed=123456,
                            plot="png")

#3.choie the ultimate clusters
clusterNum=4
cluster=result[[clusterNum]][["consensusClass"]]
write.table(cluster, file="cluster.txt", sep="\t", quote=F, col.names=F)
