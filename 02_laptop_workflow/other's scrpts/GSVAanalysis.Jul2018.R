#library(limma)
library(GSVA)
library(GSVAdata)
library(parallel)

##First read in the arguments listed at the command line
args=(commandArgs(trailingOnly=TRUE))

##args is now a list of character vectors
## First check to see if arguments are passed.
## Then cycle through each element of the list and evaluate the expressions.
if(length(args)==0){
    print("No arguments supplied.")
    stop()
}else{
    for(i in 1:length(args)){
         eval(parse(text=args[[i]]))
    }
}

#Example
#geneSetList = "/home/lolab/Downloads/h.all.v6.0.symbols.size15.gmt"

# download gene set :
# https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp#C2
# https://www.gsea-msigdb.org/gsea/msigdb/download_file.jsp?filePath=/msigdb/release/7.0/h.all.v7.0.symbols.gmt
# Download GMT Files
# gene symbols
# We want H, the CGP subset for C2, C6, and C7.
# C7 is only run when your samples have immune cells.
# So for the copied 24 samples, they are cell lines so they don't need C7
# For the download samples, they are patient biopsy, so they need C7

# for RNA-seq the input files are input.htseqCount
# if you don't have raw count matrix, you have to use the microarray mode, so input files are the same for the infile and nExpr,

# If you have the raw count, the input file is input.htseqCount and .logCPM (which is calulated, if you don't have the calculated one)

#infile = "/home/lolab/Downloads/input.htseqCount"
#nExpr = "/home/lolab/Downloads/input.logCPM"

geneSets <- getGmt(geneSetList)

# raw count from htseq-count
exprData <- read.delim(infile,sep = "\t", header = T, row.names=1)
dim(exprData)
# log-normalized exression matrix, can be cpm or fpkm.
# we need to convert the raw count to cpm, which is easier.
# conversion of raw count to cpm is using edgeR cpm
# https://www.rdocumentation.org/packages/edgeR/versions/3.14.0/topics/cpm

normExprData <- read.delim(nExpr,sep = "\t", header = T, row.names=1)
dim(normExprData)

# detect log or not logged.
# the normalized expression data with cpm must be logged with base 2.
if(transformToLog==1){
	normExprData=log2(normExprData+0.01)
}
##remove genes whose normalized expressions IQR is less than 2 fold (IQR in log2 < 1)
# remove genes that don't have big variation of expression change among whole samples.
# Apply several filters as follows:
normExprData$IQR <- apply(normExprData,1,IQR)
normExprData$max <- apply(normExprData,1,max)
normExprData<-normExprData[normExprData$IQR>=1,]
normExprData<-normExprData[normExprData$max>=minExp,]

# use normalized marix to do the subset on the raw matrix.
exprData<-exprData[rownames(normExprData),,drop=T]
dim(exprData)

##this is optional
outtable=paste(outdir,"/","EXPFilt.xls",sep="")
write.table(exprData,sep="\t",file=outtable)
#geneSets

## estimate GSVA enrichment scores
if(mArray==1){
	##for microarray data or normally distributed expression
	gsva_es <- gsva(as.matrix(exprData), geneSets, min.sz=10, max.sz=500, verbose=TRUE)
}else{
	##for RNAseq count data, Poisson distr
	# we use poisson
	gsva_es <- gsva(as.matrix(exprData), geneSets, min.sz=10, max.sz=500, kcdf="Poisson", verbose=TRUE)
}

##write the GSVA score
outfile1=paste(outprefix,"ES.xls",sep="")
write.table(gsva_es,file=outfile1,sep='\t')


