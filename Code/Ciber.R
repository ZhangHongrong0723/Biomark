#install.packages("corrplot")
#install.packages("pheatmap")


library(pheatmap)
library(corrplot)

inputFile="CIBERSORT-Results.txt"
setwd(path_way)

rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)


con=grepl("_Control", rownames(rt), ignore.case=T)
treat=grepl("_Treat", rownames(rt), ignore.case=T)
conData=rt[con,]
treatData=rt[treat,]
conNum=nrow(conData)
treatNum=nrow(treatData)
data=t(rbind(conData,treatData))


data=data[apply(data,1,sd)>0,]
Type=c(rep("Control",conNum), rep("Treat",treatNum))
names(Type)=colnames(data)
Type=as.data.frame(Type)

pdf(file="heatmap.pdf", width=7, height=4)
pheatmap(data, 
         annotation_col=Type, 
         color=colorRampPalette(c(rep("blue",3), "white", rep("red",3)))(50),
         cluster_cols=F,
         show_colnames=F,
         scale="row",
         fontsize = 7,
         fontsize_row=7,
         fontsize_col=7)
dev.off()
# col2 = colorRampPalette(c('#7F0000', 'red', '#FF7F00', 'yellow', 'white',
#                           'cyan', '#007FFF', 'blue', '#00007F'))

treatData=treatData[,apply(treatData,2,sd)>0]
pdf(file="corHeatmap.pdf", width=15, height=15,family ="sans" )
corrplot(corr=cor(treatData),
         method = "circle", 
         type="upper",
         order = "hclust", 
         tl.pos = "lt",
         insig="blank",
         sig.level = 0.1,
         pch.cex = 0.7,
         tl.col="black",
         tl.cex=1.5,
         tl.srt = 45,
         tl.offset=0.5,
         cl.cex=1.5
         # number.cex = 0.8,          
         # addCoef.col = "black",     
)
corrplot(
  corr=cor(treatData),
  add=TRUE,
  method="number",
  type="lower",
  tl.pos = "n",
  cl.pos = "n",
  diag=FALSE ,
  # col=colorRampPalette(c("#5296cc","#fdedf6","#f064af"))(50),
)
dev.off()