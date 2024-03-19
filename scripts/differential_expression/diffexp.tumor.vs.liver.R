########################################################################

library(DESeq2)
library(BiocParallel)

register(MulticoreParam(4))

########################################################################

pathExpression="../../results/expression_estimation/"
pathDifferentialExpression="../../results/differential_expression/"
pathRData="../../data_for_publication/RData/"

annot="AllTranscripts_Ensembl109_noMT_norRNA_nohaplo"

########################################################################

load(paste(pathRData, "data.sample.info.RData",sep=""))

########################################################################

read.counts=read.table(paste(pathExpression, annot, "/AllSamples_KallistoEstimatedCounts.txt",sep=""), h=T, stringsAsFactors=F)

read.counts=round(read.counts)

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

samples=c(tumor.samples$tumor_biopsyID, liver.samples$biopsyID)

files=paste(pathExpression, "/", annot,"/", samples, "/abundance.h5", sep="")
names(files)=samples

txi.kallisto <- tximport(files, type = "kallisto", tx2gene = tx2gene)

########################################################################
## prepare col data

sex=c(tumor.samples$sex, liver.samples$sex)

tissue=c(rep("Tumor", nrow(tumor.samples)), rep("Liver", nrow(liver.samples)))

colData=data.frame("Sex"=as.factor(sex), "Tissue"=as.factor(tissue))

dds=DESeqDataSetFromMatrix(read.counts, colData=colData, design = ~Sex+Tissue)

dds=DESeq(dds, test="Wald",  minReplicatesForReplace=100, parallel=T)

########################################################################

results=results(dds, contrast=c("Tissue", "Tumor", "Liver"))

results=lfcShrink(dds, coef="Tissue_Tumor_vs_Liver", res=results, type="apeglm")

results=as.data.frame(results)

results=results[order(results$padj),]

########################################################################

write.table(results, file=paste(pathDifferentialExpression, annot, "/DifferentialExpression_Tumor_vs_Liver.txt",sep=""), row.names=T, col.names=T, sep="\t", quote=F)

########################################################################
