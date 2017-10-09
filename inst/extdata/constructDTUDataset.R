local=FALSE
if(local){
setwd("/Users/koenvandenberge/PhD_Data/dtu/humanCancer/prostateCancer/")
library(tidyr)
dataMessy <- read.csv(file="/Users/koenvandenberge/PhD_Data/dtu/humanCancer/prostateCancer/kallisto_table_unnormalized_unfiltered.csv",header=TRUE)
dataMessy <- dataMessy[,c("target_id","est_counts","sample")]
dataClean <- tidyr::spread(dataMessy,key=sample,value=est_counts)
rm(dataMessy)
rownames(dataClean) <- dataClean[,"target_id"]
data <- dataClean[,-1]
prostateData=data[1:10000,]
library(biomaRt)
mart = biomaRt::useMart("ENSEMBL_MART_ENSEMBL", host='mar2016.archive.ensembl.org', dataset="hsapiens_gene_ensembl")
tx2gene = getBM(attributes=c("ensembl_transcript_id","ensembl_gene_id"), filters="ensembl_transcript_id", values=rownames(prostateData), mart=mart, bmHeader=TRUE, uniqueRows=TRUE)
colnames(tx2gene) <- c("Ensembl.Transcript.ID","Ensembl.Gene.ID")

### clean up
prostateData <- ceiling(prostateData)
## remove tx without gene match
prostateData <- prostateData[!is.na(match(rownames(prostateData),tx2gene$Ensembl.Transcript.ID)),]
## remove all zero rows
prostateData <- prostateData[!rowSums(prostateData)==0,]
## remove genes with only one tx
geneTable <- table(as.character(tx2gene$Ensembl.Gene.ID[match(rownames(prostateData),tx2gene$Ensembl.Transcript.ID)]))
genesWithOneTx <- names(geneTable)[geneTable==1]
txFromGenesWithOneTx <- tx2gene$Ensembl.Transcript.ID[tx2gene$Ensembl.Gene.ID%in%genesWithOneTx]
prostateData <- prostateData[!rownames(prostateData)%in%as.character(txFromGenesWithOneTx),]

txGeneData = as.data.frame(cbind(rownames(prostateData),as.character(tx2gene$Ensembl.Transcript.ID[match(rownames(prostateData),tx2gene$Ensembl.Transcript.ID)]),as.character(tx2gene$Ensembl.Gene.ID[match(rownames(prostateData),tx2gene$Ensembl.Transcript.ID)])))
txGeneData=txGeneData[,2:3]
colnames(txGeneData)=c("transcript","gene")
rownames(txGeneData)=txGeneData[,"transcript"]
barplot(table(table(txGeneData$gene)), main="Distribution of number of tx per gene")

#this leaves us with
length(unique(txGeneData$gene)) #nr genes
median(table(as.character(txGeneData$gene))) #median nr of tx/gene

## metadata
metaData <- read.table("sampleDataRelationship.txt",header=TRUE,sep="\t")
assays <- metaData$Assay.Name
runs <- as.character(metaData$Comment.ENA_RUN.)[seq(1,length(assays),2)]
samples=gsub(x=assays,pattern="_[1-2]",replacement="")[seq(1,length(assays),2)]
patient=factor(sapply(samples,function(x) substr(x,1,nchar(x)-1)))
condition=factor(sapply(samples,function(x) substr(x,nchar(x),nchar(x))))
prostateData <- prostateData[,match(runs,colnames(prostateData))] #same ordering as metadata
sampleData <- data.frame(condition=condition,patient=patient)
rownames(sampleData)=colnames(prostateData)


### build an expressionset.
library(Biobase)
esetProstate = ExpressionSet(assayData=as.matrix(prostateData),
			     phenoData=AnnotatedDataFrame(sampleData),
			     featureData=AnnotatedDataFrame(txGeneData))
save(esetProstate,file="~/esetProstate.RData")
devtools::use_data(esetProstate)
}
