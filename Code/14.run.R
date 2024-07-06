
inputFile="normalize.txt"      
refFile="humanREF.txt"        
setwd(path_way)     
source("geoMR14.CIBERSORT.R")      

outTab=CIBERSORT(refFile, inputFile, perm=1000)

outTab=outTab[outTab[,"P-value"]<0.05,]
outTab=as.matrix(outTab[,1:(ncol(outTab)-3)])
outTab=rbind(id=colnames(outTab),outTab)
write.table(outTab, file="CIBERSORT-Results.txt", sep="\t", quote=F, col.names=F)

