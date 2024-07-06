
library(limma)
library(ggpubr)

expFile="geneMatrix.txt"
conFile="sample1.txt"
treatFile="sample2.txt"
geneFile="intersect.txt"

setwd(path_way)

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)

con=read.table(conFile, header=FALSE, sep="\t", check.names=F)
treat=read.table(treatFile, header=FALSE, sep="\t", check.names=F)
conData=data[,as.vector(con[,1])]
treatData=data[,as.vector(treat[,1])]
data=cbind(conData,treatData)
conNum=ncol(conData)
treatNum=ncol(treatData)

Type=c(rep("con",conNum),rep("treat",treatNum))
outData=rbind(id=paste0(colnames(data),"_",Type),data)
# write.table(outData,file="result.txt",sep="\t",check.names=FALSE)
geneRT=read.table(geneFile, header=FALSE, sep="\t", check.names=F)
data=data[as.vector(geneRT[,1]),,drop=FALSE]

  
Type=c(rep("Con",conNum),rep("Treat",treatNum))
mycomparised=list()
mycomparised[[1]]=levels(factor(Type))
newGeneLists=c()
outTab=data.frame()
for (i in row.names(data)) {
  rt1=data.frame(expression=data[i,],Type=Type)
  boxplot=ggboxplot(rt1, x="Type", y="gene", color="Type",
                    xlab="",
                    ylab=paste(i, " expression"),
                    legend.title="Type",
                    palette = c("blue","red"),
                    add = "jitter")+ 
    stat_compare_means(comparisons=mycomparised)
  pdf(file=paste0("boxplot.",i,".diff.pdf"), width=5, height=4.5)
  print(boxplot)
  dev.off()
  
}