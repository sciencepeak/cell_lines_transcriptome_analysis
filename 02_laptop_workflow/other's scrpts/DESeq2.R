library('DESeq2')
library("ggplot2")

#load table
countsTable<-read.delim("AR12202018.htseq.txt",sep="\t",header=T,row.names=1)
colnames(countsTable)
#reorder the column by name
t=countsTable[order(colnames(countsTable))]

#Design
Treatment<-factor(c(rep("DMSO.model1",2),rep("DMSO.model2",3),rep("Vem.model1",4),rep("Vem.model2",3)),
                  levels=c("DMSO.model1","Vem.model1","DMSO.model2","Vem.model2"))
des<-formula(~Treatment)

#Make a table assoc. the sample name and treatment
myNames<-colnames(t)
colDataNames<-data.frame(row.names=myNames, Treatment=Treatment)
ddsHTSeq<-DESeqDataSetFromMatrix(t, colData=colDataNames, design=des, ignoreRank = FALSE)

dds<-DESeq(ddsHTSeq,betaPrior=FALSE)
resultsNames(dds)

normCounts<-as.data.frame(counts(dds,normalized=TRUE))
write.csv(normCounts,"AR12202018.DESeq_outputNormalised.csv")

#PCA 
rld <- rlogTransformation(dds, blind=TRUE)
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
PCAdata<-plotPCA(rld, intgroup=c("Treatment"),ntop = 1000,returnData=TRUE)
percentVar <- round(100 * attr(PCAdata, "percentVar"))
myplot1<-qplot(PC1, PC2, color=Treatment, data=PCAdata) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
ggtitle("PCA All Conditions") +
scale_size_manual(values=c(3,5,7,9))+
geom_point(size=6)+
guides(colour = guide_legend(override.aes = list(size=6))) +
guides(shape = guide_legend(override.aes = list(size=6))) +
theme(plot.title = element_text(colour="black", size = 16)) + 
theme(axis.title = element_text(colour="black", size = 14)) + 
theme(legend.title = element_text(colour="black", size = 16)) + 
theme(legend.text = element_text(colour="black", size = 16)) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
myplot1




#DE after Removal Control 4

res<-results(dds,contrast=list("Treatment_MEKiLast_vs_MEKiCont"),cooksCutoff=FALSE,independentFiltering=FALSE)
write.csv(res,"DESeq_Treatment_MEKiCont_vs_MEKiLast.Dec.2016.csv")

res<-results(dds,contrast=list("Treatment_BRAFiLast_vs_MEKiCont"),cooksCutoff=FALSE,independentFiltering=FALSE)
write.csv(res,"DESeq_Treatment_MEKiCont_vs_BRAFiLast.Dec.2016.csv")

res<-results(dds,contrast=list("Treatment_MEKiWD_vs_MEKiCont"),cooksCutoff=FALSE,independentFiltering=FALSE)
write.csv(res,"DESeq_Treatment_MEKiWD_vs_MEKiCont.Dec.2016.csv")

res<-results(dds,contrast=list("Treatment_BRAFi_vs_MEKiCont"),cooksCutoff=FALSE,independentFiltering=FALSE)
write.csv(res,"DESeq_Treatment_BRAFi_vs_MEKiCont.Dec.2016.csv")

Treatment<-factor(c(rep("MEKiLast",3),rep("BRAFiLast",3),rep("MEKiCont",4),rep("MEKiWD",4),rep("BRAFi",4)),levels=c("MEKiLast","BRAFiLast","BRAFi","MEKiCont","MEKiWD"))
des<-formula(~Treatment)

#Data
countsTable<-read.delim("RL11092016_RNAseq_MTGR2_5mg_BRAFi.HtseqCount",sep="\t",header=T,row.names=1)
myNames<-colnames(countsTable)[1:18]
colDataNames<-data.frame(row.names=myNames, Treatment=Treatment)
ddsHTSeq<-DESeqDataSetFromMatrix(countsTable, colData=colDataNames, design=des, ignoreRank = FALSE)
dds<-DESeq(ddsHTSeq,betaPrior=FALSE)
resultsNames(dds)
normCounts<-as.data.frame(counts(dds,normalized=TRUE))
write.csv(normCounts,"RL11092016_RNAseq_MTGR2_5mg_BRAFi.DESeq_outputNormalised.csv")



