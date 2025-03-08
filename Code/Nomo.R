
library(rms)
library(rmda)

inputFile="normalize.txt"    
geneFile="hubGenes.txt"       
setwd(path_way)     

data=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
row.names(data)=gsub("-", "_", row.names(data))

RT=read.table(geneFile, header=T, sep="\t", check.names=F)
data=data[as.vector(geneRT[,1]),]

row.names(data)){
	data[i,]=ifelse(as.numeric(data[i,])>median(as.numeric(data[i,])), "High", "Low")
}

group=gsub("(.*)\\_(.*)", "\\2", row.names(data))
rt=cbind(as.data.frame(data), Type=group)


options(datadist="ddist")

nomo=nomogram(lrmModel, fun=plogis,
	fun.at=c(0.001,0.1,0.5,0.9,0.99),
	lp=F, funlabel="Risk of Disease")

plot(8,family="sans"omo)
dev.off()


pdf("Calibration.pdf", width=5.5, height=5.5)
plot(c,family="sans"ali,
	xlab="Predicted probability",
	ylab="Actual probability", sub=F)
dev.off()
