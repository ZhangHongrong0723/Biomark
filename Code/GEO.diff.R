
library(limma)
library(pheatmap)
library(ggplot2)

inputFile="GSE7014.txt"   
conFile="GSE7014_s1.txt"             
treatFile="GSE7014_s2.txt"           
logFCfilter=1                  
adj.P.Val.Filter=0.01         
setwd(path_way)     

rt=read.table(inputFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
rt=avereps(data)

qx=as.numeric(quantile(rt, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC=( (qx[5]>100) || ( (qx[6]-qx[1])>50 && qx[2]>0) )
if(LogC){
	rt[rt<0]=0
	rt=log2(rt+1)}
data=normalizeBetweenArrays(rt)

sample1=read.table(conFile, header=F, sep="\t", check.names=F)
sample2=read.table(treatFile, header=F, sep="\t", check.names=F)
sampleName1=gsub("^ | $", "", as.vector(sample1[,1]))
sampleName2=gsub("^ | $", "", as.vector(sample2[,1]))
conData=data[,sampleName1]
treatData=data[,sampleName2]
data=cbind(conData,treatData)
conNum=ncol(conData)
treatNum=ncol(treatData)

Type=c(rep("con",conNum), rep("treat",treatNum))
design <- model.matrix(~0+factor(Type))
colnames(design) <- c("con","treat")
fit <- lmFit(data,design)
cont.matrix<-makeContrasts(treat-con,levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

allDiff=topTable(fit2, adjust='BY', number=200000)
allDiffOut=rbind(id=colnames(allDiff),allDiff)
write.table(allDiffOut, file="all.txt", sep="\t", quote=F, col.names=F)

Type=c(rep("Control",conNum),rep("Treat",treatNum))
outData=rbind(id=paste0(colnames(data),"_",Type),data)
write.table(outData, file="normalize.txt", sep="\t", quote=F, col.names=F)

# diffSig=allDiff[with(allDiff, (abs(logFC)>logFCfilter & adj.P.Val < adj.P.Val.Filter )),]
diffSig=allDiff[with(allDiff, (logFC>1 &  P.Value<0.01)), ]
diffSigOut=rbind(id=colnames(diffSig),diffSig)
write.table(diffSigOut, file="diff.txt", sep="\t", quote=F, col.names=F)


geneNum=50  
diffUp=diffSig[diffSig$logFC>0,]
diffDown=diffSig[diffSig$logFC<0,]
geneUp=row.names(diffUp)
geneDown=row.names(diffDown)
if(nrow(diffUp)>geneNum){geneUp=row.names(diffUp)[1:geneNum]}
if(nrow(diffDown)>geneNum){geneDown=row.names(diffDown)[1:geneNum]}
hmExp=data[c(geneUp,geneDown),]

Type=c(rep("Control",conNum),rep("Treat",treatNum))
names(Type)=colnames(data)
Type=as.data.frame(Type)

pdf(file="heatmap.pdf", width=8, height=6.5)
pheatmap(hmExp, 
         annotation_col =Type, 
         color = colorRampPalette(c("blue", "white", "red"))(50),
         cluster_cols =F,
         show_colnames = F,
         scale="row",
         fontsize = 7,
         fontsize_row=5,
         fontsize_col=7)
dev.off()


allDiff$logFC[allDiff$logFC>20]=20
allDiff$logFC[allDiff$logFC< -20]=-20
# Significant=ifelse((allDiff$adj.P.Val<adj.P.Val.Filter & abs(allDiff$logFC)>logFCfilter), ifelse(allDiff$logFC>logFCfilter,"Up","Down"), "Not")
# Significant=ifelse((allDiff$adj.P.Val<0.001 & abs(allDiff$logFC)>logFCfilter & allDiff$P.Value<0.001), ifelse(allDiff$logFC>logFCfilter,"Up","Down"), "Not")
p = ggplot(allDiff,aes(logFC, -log10(P.Value)))+
    geom_hline(yintercept = -log10(0.05),linetype="dashed",color= "#999999")+
    geom_vline(xintercept = c(-1.2,1.2),linetype="dashed",color= "#999999")+
    geom_point(aes(size=-log10(P.Value),color=-log10(P.Value)))+
    labs(title="   ")+
    scale_color_gradientn(values=seq(0,1,0.2),
                        colors=c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25"))+
    scale_size_continuous(range=c(1,3))+
    theme(
      axis.text = element_text(face="bold", hjust=0.5,size=10),
      panel.border = element_rect(fill = NA, colour = "black",linewidth = 1),
      axis.title = element_text(face="bold", hjust=0.5,size=10),
      legend.text = element_text(face="bold", hjust=0.5,size=10),
      legend.title = element_text(face="bold", hjust=0.5,size=10)
    )+
    theme(panel.grid = element_blank())

# Fon <- 'sans'
# p = ggplot(allDiff, aes(logFC, -log10(adj.P.Val)),family=Fon,size=20)+
#     geom_point(aes(col=Significant))+
#     scale_color_manual(values=c("green", "black", "red"))+
#     labs(title = "   ")+
#     theme(panel.background = element_rect(fill = "white", colour = NA),
#           panel.border = element_rect(fill = NA, colour = "black",linewidth = 2),
#           plot.title = element_text(family = Fon,size=40, hjust=0.5, face = "bold"),
#           axis.text = element_text(family = Fon,face="bold", hjust=0.5,size=10) )

# p=p+theme_bw()

pdf(file="vol.pdf", width=5.5, height=5)
print(p)
dev.off()


