library(edgeR)
library(limma)
library(gplots)

infile <- "htseq_count_expression_matrix.csv.bz2"
dat <- read.csv(infile,header=T,row.names=1)

##testing to show the raw count of MITF
dat["MITF",]

#normalization into CPM
dge <- DGEList(counts=dat)
norm.factor <- calcNormFactors(dge)
dat.norm <- cpm(norm.factor, normalized.lib.sizes=TRUE, prior.count=2, log=TRUE)
#outfile <- "htseq_count_expression_matrix.CPM.txt"
#write.table(dat.norm,file=outfile,sep="\t",col.names=TRUE,row.names=TRUE)

colnames(dat.norm)

##concatenate the cpm expression tables
##remove batch effect
library(readxl)
library(stringr)
batchTable = read_xlsx("samples_cell_lines_batch_both_SampleGroup.xlsx")
batchTable$File = str_replace_all(batchTable$File,"-",".")
rownames(batchTable) = batchTable$File
batchTable = batchTable[colnames(dat.norm),]

dat.batch.adj <- removeBatchEffect(dat.norm,batch=batchTable$Batch)
write.table(dat.batch.adj,file="htseq_count_expression_matrix.CPM.batchAdj.txt",sep="\t",col.names=TRUE,row.names=TRUE)
dat.norm = as.data.frame(dat.batch.adj)


## Find pseudo genes and genes that are not confirmed
genes = rownames(dat.norm)
genesFiltered=grep("^[A-Z][A-Z1-9]+[.][1-9]+", genes, perl=TRUE, value=TRUE, invert=TRUE)
genesFiltered=grep("\\.", genesFiltered, perl=TRUE, value=TRUE, invert=TRUE)
genesFiltered=grep("^MT-*[A-Z1-9]+", genesFiltered, perl=TRUE, value=TRUE, invert=TRUE)
genesFiltered=grep("hsa-mir", genesFiltered, perl=TRUE, value=TRUE, invert=TRUE)
genesFiltered=grep("^MIR\\d+", genesFiltered, perl=TRUE, value=TRUE, invert=TRUE)
genesFiltered=grep("^RNA\\d+[A-Z1-9]+", genesFiltered, perl=TRUE, value=TRUE, invert=TRUE)
genesFiltered=grep("^RNU\\d+", genesFiltered, perl=TRUE, value=TRUE, invert=TRUE)
genesFiltered=grep("^RNVU\\d+", genesFiltered, perl=TRUE, value=TRUE, invert=TRUE)
genesFiltered=grep("^RNY\\d+", genesFiltered, perl=TRUE, value=TRUE, invert=TRUE)
genesFiltered=grep("^HNRN", genesFiltered, perl=TRUE, value=TRUE, invert=TRUE)
genesFiltered=grep("^SNOR[A-Z]\\d+[A-Z1-9]*", genesFiltered, perl=TRUE, value=TRUE, invert=TRUE)
genesFiltered=grep("^LOC\\d+[A-Z1-9]+", genesFiltered, perl=TRUE, value=TRUE, invert=TRUE)
genesFiltered=grep("^LINC\\d+[A-Z1-9]+", genesFiltered, perl=TRUE, value=TRUE, invert=TRUE)
genesFiltered=grep("^SCARNA\\d+[A-Z1-9]*", genesFiltered, perl=TRUE, value=TRUE, invert=TRUE)
genesFiltered=grep("[A-Z1-9]+-[A-Z1-9][A-Z1-9]+", genesFiltered, perl=TRUE, value=TRUE, invert=TRUE)

## remove genes that are pseudo genes/unconfirmed
dat.norm=as.data.frame(dat.norm[genesFiltered,,drop=T])

## save your intermediate
saveRDS(dat.norm,file="GenesFiltered.rds")
dat.norm = readRDS("GenesFiltered.rds")

## Identify genes with high variability
## compute IQR, max value
t=dat.norm
t$IQR=apply(dat.norm, 1, IQR)
t$max=apply(dat.norm, 1, max)

#sort by IQR
t=t[order(t$IQR,decreasing=T),]

## remove gene with log2 CPM < 0, and also choose those with IQR >= 1
t=t[t$max>=0,]
t=t[t$IQR>=1,]

columnLen = length(colnames(t))
c1 = c(-(columnLen-1),-columnLen)

##remove IQR and max column, choose the top 5000 genes based on IQR
t5000=t[1:5000,c(-(columnLen-1),-columnLen)]

##remove M238-Ctrl (outlier)
t5000=t5000[,c(-15)]

## the PCA computation, use scaling
pca.res<-prcomp(t(t5000),scale.=T,retx=T)
pc.var<-pca.res$sdev^2

#percentage of variance explained 
pc.per<-round(pc.var/sum(pc.var)*100, 1)
x<-as.data.frame(pca.res$x)

plot.default(x$PC1,x$PC2,type="p",pch=c(21),
             bg="red",
             col="white",
             lwd=2,cex=2.5,xaxs="r",yaxs="r",bty="l",
             cex.axis=1,cex.main=1,cex.lab=1,font.lab=2,font.axis=2)


##Create a copy of x with additional column "File"
x1 = x
x1$File = rownames(x1)
##### up to here###
batchTable = batchTable[-c(15),]

### noted batch effect in original data
x2 = merge(x1,batchTable,by="File")

plot.default(x2$PC1,x2$PC2,type="p",pch=c(21),
             bg=c("red","black")[x2$Batch],
             col="white",
             lwd=2,cex=2.5,xaxs="r",yaxs="r",bty="l",
             cex.axis=1,cex.main=1,cex.lab=1,font.lab=2,font.axis=2)

library(rgl)
plot3d(x2$PC1,x2$PC2,x2$PC3,type="s",col=c("red","black")[x2$Batch],size=2)
grid3d("x")
grid3d("y+")
grid3d("z")
rgl.snapshot("Combined.BatchAdj.png")


#t2 <- as.data.frame(cbind(x[,1:5],MouseModel,Treatment,Batch))

t2 <- as.data.frame(x[,1:5])
t2$MouseModel = MouseModel
t2$Treatment = Treatment
t2$Batch = Batch


xlab1=paste("PC1 ",pc.per[[1]],"%",sep="")
ylab1=paste("PC2 ",pc.per[[2]],"%",sep="")
zlab1=paste("PC3 ",pc.per[[3]],"%",sep="")

mainLabel="PCA analysis"
xlim1=c(-60,70)
#ylim1=c(-20,60)

##write the plot into TIFF
tiff(filename="PCA.Batch.tiff")
plot.default(t2$PC1,t2$PC2,type="p",pch=c(21,0,23)[factor(t2$Treatment)],
                            bg=colGrp[factor(t2$MouseModel)],
                            col="white",
                            lwd=2,cex=2.5,xlab=xlab1,ylab=ylab1,xaxs="r",yaxs="r",bty="l",
                            xlim=xlim1,cex.axis=1,cex.main=1,cex.lab=1,font.lab=2,font.axis=2)

plot.default(t2$PC1,t2$PC2)

box(lwd=3,bty="l")
axis(1,lwd=3,labels=F)
axis(2,lwd=3,labels=F)

legend(x="topright",levels(factor(Treatment)),pch=c(21:23),pt.bg=rep("black",3),pt.cex=2.5,pt.lwd=2,col="black",cex=1.2,text.font=2,bty="n")
legend(x="bottomright",levels(factor(MouseModel)),pch=20,pt.bg=colGrp,pt.cex=2.5,pt.lwd=2,col=colGrp,cex=1,text.font=2,bty="n")

dev.off()


