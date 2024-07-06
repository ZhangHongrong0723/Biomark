#install.packages("VennDiagram")
library(VennDiagram)             
outFile="intersectGenes.txt"      
outPic="venn.pdf"                  
setwd(path_way)   
files=dir()                        
files=grep("txt$",files,value=T)   
geneList=list()


for(i in 1:length(files)){
    inputFile=files[i]
	if(inputFile==outFile){next}
    rt=read.table(inputFile,header=F)        
    geneNames=as.vector(rt[,1])            
    geneName=toupper(geneNames)
    geneNames=gsub("^ | $","",geneNames)     
    uniqGene=unique(geneNames)              
    header=unlist(strsplit(inputFile,"\\.|\\-"))
    geneList[[header[1]]]=uniqGene
    uniqLength=length(uniqGene)
    print(paste(header[1],uniqLength,sep=" "))
}


color=c("#87CEFA", "#ff6961","#98fb98")#
venn.plot=venn.diagram(geneList,filename=NULL,fill=color,lty = 'blank',sub.cex = 2)
pdf(file=outPic, width=5, height=5,family="sans")
grid.draw(venn.plot)
dev.off()


intersectGenes=Reduce(intersect,geneList)
write.table(file=outFile,intersectGenes,sep="\t",quote=F,col.names=F,row.names=F)
