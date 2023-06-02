#install.packages("survival")
#install.packages("caret")
#install.packages("glmnet")
#install.packages("survminer")
#install.packages("timeROC")



library(survival)
library(caret)
library(glmnet)
library(survminer)
library(timeROC)
setwd("C:\\Users\\lenovo\\Desktop\\lasso analysis") 
#1. read file     
rt=read.table("uniSigExp.txt", header=T, sep="\t", check.names=F, row.names=1)     
rt$futime[rt$futime<=0]=0.003

#2. set the train and test group 
for(i in 1:1000){
	
	inTrain<-createDataPartition(y=rt[,3], p=0.5, list=F)
	train<-rt[inTrain,]
	test<-rt[-inTrain,]
	trainOut=cbind(id=row.names(train), train)
	testOut=cbind(id=row.names(test), test)

	#3. lasso regression 
	x=as.matrix(train[,c(3:ncol(train))])
	y=data.matrix(Surv(train$futime,train$fustat))
	fit <- glmnet(x, y, family = "cox", maxit = 1000)
	cvfit <- cv.glmnet(x, y, family="cox", maxit = 1000)
	coef <- coef(fit, s = cvfit$lambda.min)
	index <- which(coef != 0)
	actCoef <- coef[index]
	lassoGene=row.names(coef)[index]
	geneCoef=cbind(Gene=lassoGene, Coef=actCoef)
	if(nrow(geneCoef)<2){next}
	
	#4. output the risk score about train group 
	trainFinalGeneExp=train[,lassoGene]
	myFun=function(x){crossprod(as.numeric(x), actCoef)}
	trainScore=apply(trainFinalGeneExp, 1, myFun)
	trainScore=predict(cvfit, newx=as.matrix(train[,c(3:ncol(train))]), s="lambda.min", type="response")
	outCol=c("futime", "fustat", lassoGene)
	risk=as.vector(ifelse(trainScore>median(trainScore), "high", "low"))
	train=cbind(train[,outCol], riskScore=as.vector(trainScore), risk)
	trainRiskOut=cbind(id=rownames(train), train)

	#5. output the risk score about test group 
	testFinalGeneExp=test[,lassoGene]
	testScore=apply(testFinalGeneExp, 1, myFun)
	testScore=predict(cvfit, newx=as.matrix(test[,c(3:ncol(test))]), s="lambda.min", type="response")
	outCol=c("futime", "fustat", lassoGene)
	risk=as.vector(ifelse(testScore>median(trainScore), "high", "low"))
	test=cbind(test[,outCol], riskScore=as.vector(testScore), risk)
	testRiskOut=cbind(id=rownames(test), test)
	
	#6. survival analysis	
	diff=survdiff(Surv(futime, fustat) ~risk, data=train)
	pValue=1-pchisq(diff$chisq, df=1)
	diffTest=survdiff(Surv(futime, fustat) ~risk, data=test)
	pValueTest=1-pchisq(diffTest$chisq, df=1)

	#7. paint the ROC
	predictTime=3    
	roc=timeROC(T=train$futime, delta=train$fustat,
	            marker=trainScore, cause=1,
	            weighting='aalen',
	            times=c(predictTime), ROC=TRUE)
	rocTest=timeROC(T=test$futime, delta=test$fustat,
	            marker=testScore, cause=1,
	            weighting='aalen',
	            times=c(predictTime), ROC=TRUE)	
	
	if((pValue<0.01) & (roc$AUC[2]>0.75) & (pValueTest<0.05) & (rocTest$AUC[2]>0.7)){
		#8.save results about grouping 
		write.table(trainOut,file="train.data.txt",sep="\t",quote=F,row.names=F)
		write.table(testOut,file="test.data.txt",sep="\t",quote=F,row.names=F)
		
		pdf("lambda.pdf")
		plot(fit, xvar = "lambda", label = TRUE)
		dev.off()
		pdf("cvfit.pdf")
		plot(cvfit)
		abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)),lty="dashed")
		dev.off()
	   
		write.table(geneCoef, file="geneCoef.txt", sep="\t", quote=F, row.names=F)
		write.table(trainRiskOut,file="trainRisk.txt",sep="\t",quote=F,row.names=F)
		write.table(testRiskOut,file="testRisk.txt",sep="\t",quote=F,row.names=F)
		
		allRiskOut=rbind(trainRiskOut, testRiskOut)
		write.table(allRiskOut,file="allRisk.txt",sep="\t",quote=F,row.names=F)
		break
	}
}


