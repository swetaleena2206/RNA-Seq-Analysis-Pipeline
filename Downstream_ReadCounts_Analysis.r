#####...IMPORTING DATA AND LOADING LIBRARIES...#####
#Loading the required libraries (some of them are from Bioconductor package)
library(edgeR)
library(limma)
library(Glimma)
library(gplots)
library(org.Hs.eg.db)
library(RColorBrewer)
library(dplyr)
library(topGO)
library(pathview)
library(AnnotationDbi)

#Reading and viewing the ReadCounts file
countdata <- read.delim("SRP002628_ReadCounts.txt", sep=',',stringsAsFactors = FALSE)
head(countdata)

#####...DATA PREPROCESSING AND CLEANING...#####
#making sure the transcripts are all added up
countdata <- aggregate(countdata[-1], countdata[1], mean)

#replacing any transcript IDs 
countdata$GeneID <- gsub("\\..*","",countdata$Gene_ID) 

#converting the processed GeneID columns to row.names
rownames(countdata) <- countdata[,1]
head(countdata)		

#annotating the REFSEQ IDs into Symbols and Genenames
annotations <-AnnotationDbi::select(org.Hs.eg.db, keys = row.names(countdata),column = c("SYMBOL","GENENAME"), keytype = "REFSEQ", multiVals ="first")
head(annotations)

#merging the annotated data into the countdata
countdata <- merge(countdata, annotations, by.x=0, by.y="REFSEQ")
head(countdata)										

#subsetting and extracting only what we require at the moment
countdata2<-countdata[,c(12,3:10)]
head(countdata2)

#further cleaning and processing of data
countdata2 <- aggregate(countdata2[-1], countdata2[1], mean)
rownames(countdata2) <- countdata2[,1]
head(countdata2)

countdata2$SYMBOL <- NULL
head(countdata2)

#Creating Groups for analysis
group <- c(rep("Tumour",4),rep("Control",4))

#####...FILTERING TO REMOVE LOWLY EXPRESSED GENES...##### 
#Calculating the counts per million value for each gene in each sample 
myCPM <- cpm(countdata2) #CPM
head(myCPM)

#Plotting a graph of CPM vs Counts to decide a threshold 
#We will look at the first sample
#There is a corresponding value of CPM = 0.65 intersecting the Counts axis at 10. 
plot(myCPM[,1],countdata2[,1], xlab="CPM", ylab="Counts", ylim=c(0,50),xlim=c(0,3) , main="CPM vs Counts",
pch=1, cex.main=1.5, frame.plot=FALSE , col="blue")
abline(v=0.65, h=10,col=c("blue", "red"))
 

#setting a CPM threshold =0.65
thresh <- myCPM > 0.65  
#This produces a logical matrix with TRUEs and FALSEs
head(thresh)

#Summary of how many TRUEs there are in each row
#There are 16591 genes that have TRUEs in all 8 samples.
table(rowSums(thresh))

#keeping those genes which have at least 2 TRUES in each row
keep <- rowSums(thresh) >= 2

#Subset the rows of countdata to keep the more highly expressed genes
countdata2.keep <- countdata2[keep,]
head(countdata2.keep)

#To check the dimensions of countdata2.keep
dim(countdata2.keep)


#####...CREATING A DGE OBJECT AFTER FIRST ROUND OF NORMALIZATION...#####
dgeObj <- DGEList(counts = countdata2.keep, group = group)

#have a look at DGE object
dgeObj

#To check what slots are stored in dgeObj
names(dgeObj)

#Library size information is stored in the samples slot
dgeObj$samples

#####...QUALITY CONTROL...#####
#plot the library sizes as a barplot to see whether there are any major discrepancies between the samples
#The names argument tells the barplot to use the sample names on the x-axis
#The las argument rotates the axis names
barplot(dgeObj$samples$lib.size,names=colnames(dgeObj),las=2)
title("Barplot of library sizes")

#Boxplot before normalization
#Get log2 counts per million
logcounts <- cpm(dgeObj,log=TRUE)
col.cell <- c("red","red","red","red","green","green","green","green")

#Check distributions of samples using boxplots
boxplot(logcounts, xlab="", ylab="Log2 counts per million",col=col.cell,las=2)
#adding a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts),col="blue")
title("Boxplots of logCPMs (unnormalised)")

#Plotting Heatmap using highest variable genes
logcounts <- cpm(dgeObj,log=TRUE)

#We estimate the variance for each row in the logcounts matrix 
var_genes <- apply(logcounts, 1, var)
head(var_genes)

#Get the gene names for the top 500 most variable genes
select_var <- names(sort(var_genes, decreasing=TRUE))[1:500] 
head(select_var)

#Subset logcounts matrix
highly_variable_lcpm <- logcounts[select_var,]
dim(highly_variable_lcpm)
head(highly_variable_lcpm)

#Get some nicer colours
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)

#Plot the heatmap 
heatmap.2(highly_variable_lcpm,col=rev(morecols(50)), trace="column",main="Top 100 most variable genes across samples", ColSideColors=col.cell,
scale="row")

#Save the heatmap
png(file="High_var_genes.heatmap.png")
heatmap.2(highly_variable_lcpm,col=rev(morecols(50)),trace="none",main="Top 100 most variable genes\nacross samples",ColSideColors=col.cell ,
scale="row")
dev.off()

#####...NORMALIZATION FOR COMPOSITIONAL BIAS...#####
#Apply normalisation to DGEList object
dgeObj <- calcNormFactors(dgeObj,method="upperquartile")

#To take a look at the normalization factors for the samples
summary(dgeObj)
head(dgeObj$samples)

#Controling for Dispersion (Biological coefficient of variation)
dgeObj <- estimateCommonDisp(dgeObj, verbose=TRUE)
dgeObj <- estimateGLMTrendedDisp(dgeObj,verbose=TRUE)
dgeObj <- estimateTagwiseDisp(dgeObj,verbose=TRUE)
#the BCV is used to either reduce or increase the dispersion value of a gene.

#####...BOX PLOT AFTER NORMALIZATION...#####
normlogcounts <- cpm(dgeObj,log=TRUE)
boxplot(normlogcounts, xlab="", ylab="Log2 counts per million",col=col.cell,las=2)
abline(h=median(normlogcounts),col="blue")
title("Boxplots of logCPMs (normalised)")

#Putting the box-plots side-by-side for comparison
par(mfrow=c(1,2))
boxplot(logcounts, xlab="", ylab="Log2 counts per million",col=col.cell,las=2,main="Unnormalised logCPM")

#Let's add a blue horizontal line that corresponds to the median logCPM 
abline(h=median(logcounts),col="blue")
boxplot(normlogcounts, xlab="", ylab="Log2 counts per million" ,col=col.cell,las=2, main="Normalised logCPM") 
abline(h=median(normlogcounts),col="blue")
 
#####...DIFFERENTIAL GENE EXPRESSION...#####
dgeanalysistest <- exactTest(dgeObj)
etp<-topTags(dgeanalysistest,n=nrow(dgeanalysistest$table))
head(etp$table)

#Multiple Correction adjusted p-value
dgeanalysistest$PValue_fdr <- p.adjust(method="fdr",p=dgeanalysistest$PValue) 
dgetable <- dgeanalysistest$table
dgetable$PValue_fdr <- p.adjust(method="fdr",p=dgetable$PValue)
head(dgetable)

#5% False Discovery rate 
table(dgetable$PValue_fdr<0.05) 
isde<-decideTestsDGE(dgeanalysistest,p.value=0.05) 
summary(isde)
#extracting only the DEGs
de_results <- as.data.frame(dgetable)
fdr_cutoff <- 0.2

de_results$de <- dgetable$PValue_fdr < fdr_cutoff 
head(de_results)
write.table(de_results, file = "de-results.txt", quote = FALSE, sep = "\t", col.names = NA) #this will be used for gene annotation

plotMD(dgeanalysistest, status=isde, values=c(1,-1), col=c("red","blue"), legend="topright")

#Annotate the expression data and saving
ann <- try(suppressWarnings(AnnotationDbi::select(org.Hs.eg.db,keys=rownames (dgetable),columns=c("ENTREZID","GENENAME"), keytype='SYMBOL')))
head(ann)

#merging tables
diffexprwithannotation<-merge(etp$table,ann,by.x="row.names",by.y="SYMBOL") 
write.table(diffexprwithannotation, file="DifferentiallyExpressedGenes.csv" , sep='\t',quote=FALSE)

#Creating Heatmap of DEGs
getgeneorder <- order(dgetable$PValue_fdr)
head(getgeneorder)
#top 30 DEGs will be plotted 
newlogCPM <- normlogcounts[getgeneorder[1:30],]
head(newlogCPM)

coolmap(newlogCPM, margins=c(7,7), lhei=c(1,6), lwid=c(1,3))