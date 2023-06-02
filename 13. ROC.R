#install.packages("survival")
#install.packages("survminer")
#install.packages("timeROC")


library(survival)
library(survminer)
library(timeROC)
setwd("C:\\Users\\lenovo\\Desktop\\26ROC")      


bioROC=function(inputFile=null, rocFile=null){
	predictTime=3   
	
	rt=read.table(inputFile, header=T, sep="\t")
	
	ROC_rt=timeROC(T=rt$futime, delta=rt$fustat,
	               marker=rt$riskScore, cause=1,
	               weighting='aalen',
	               times=c(predictTime), ROC=TRUE)
	pdf(file=rocFile, width=5, height=5)
	plot(ROC_rt, time=predictTime, col='red', title=FALSE, lwd=2)
	legend('bottomright', cex=1.3,
           paste0('AUC=',sprintf("%.03f",ROC_rt$AUC[2])),
	       col="white", lwd=1, bty = 'n')
	dev.off()
}

bioROC(inputFile="trainRisk.txt",rocFile="train.ROC.pdf")
bioROC(inputFile="testRisk.txt",rocFile="test.ROC.pdf")


