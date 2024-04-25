########################################################################

library(DESeq2)
library(BiocParallel)
library(tximport)

register(MulticoreParam(4))

########################################################################

pathExpression="../../results/expression_estimation/"
pathDifferentialExpression="../../results/differential_expression/"
pathRData="../../data_for_publication/RData/"

annot="AllTranscripts_Ensembl109_noMT_norRNA_nohaplo"

########################################################################

load(paste(pathRData, "data.sample.info.MERiC.RData",sep=""))

########################################################################

## get read counts for each gene

read.counts=read.table(paste(pathExpression, annot, "/AllSamples_KallistoEstimatedCounts.txt",sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")

total.estimated.counts=apply(read.counts,2,sum)

tumor.samples$total.estimated.counts=total.estimated.counts[tumor.samples$tumor_biopsyID]

## order samples by read count and Edmondson degree
tumor.samples=tumor.samples[order(tumor.samples$total.estimated.counts, decreasing=T),]

tumor.samples=tumor.samples[order(tumor.samples$edmondson, decreasing=T),]

## select one sample per patient

dupli=which(duplicated(tumor.samples$Patient_ID))

if(length(dupli)>0){
    tumor.samples=tumor.samples[-dupli,]
}

########################################################################

## get tx2gene info

samples=tumor.samples$tumor_biopsyID

sinfo=read.table(paste(pathExpression, annot, "/",samples[1],"/abundance.tsv",sep=""), h=T, stringsAsFactors=F)
geneid=unlist(lapply(sinfo$target_id, function(x) unlist(strsplit(x, split=":"))[1]))

tx2gene=data.frame("txid"=sinfo$target_id, "geneid"=geneid)

########################################################################

## assemble kallisto counts

files=paste(pathExpression, "/", annot,"/", samples, "/abundance.h5", sep="")
names(files)=samples

txi.kallisto <- tximport(files, type = "kallisto", tx2gene = tx2gene)

########################################################################

for(disease in c("
