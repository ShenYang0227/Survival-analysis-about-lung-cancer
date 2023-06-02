#install.packages("survival")
#install.packages("survminer")



library(survival)
library(survminer)
setwd("C:\\Users\\lenovo\\Desktop\\25 survival")       #设置工作目录


bioSurvival=function(inputFile=null,outFile=null){
	
	rt=read.table(inputFile, header=T, sep="\t")
	
	diff=survdiff(Surv(futime, fustat) ~risk,data = rt)
	pValue=1-pchisq(diff$chisq,df=1)
	if(pValue<0.001){
		pValue="p<0.001"
	}else{
		pValue=paste0("p=",sprintf("%.03f",pValue))
	}
	fit <- survfit(Surv(futime, fustat) ~ risk, data = rt)
		
	
	surPlot=ggsurvplot(fit,pval = T,
           conf.int = T,
           conf.int.style="ribbon",
           surv.median.line = "hv",
           conf.int.alpha=0.1,
           palette = "jco",
           risk.table = T)

	pdf(file=outFile,onefile = FALSE,width = 6.5,height =5.5)
	print(surPlot)
	dev.off()
}
bioSurvival(inputFile="trainRisk.txt", outFile="trainSurv.pdf")
bioSurvival(inputFile="testRisk.txt", outFile="testSurv.pdf")


