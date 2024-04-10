########################################################################

library(DESeq2)
library(BiocParallel)
library(tximport)

register(MulticoreParam(4))

########################################################################

pathExpression="../../results/expression_estimation_TCGA/"
pathDifferentialExpression="../../results/differential_expression_TCGA/"
pathRData="../../data_for_publication/RData/"

annot="AllTranscripts_Ensembl109_noMT_norRNA_nohaplo"

########################################################################

load(paste(pathRData, "data.sample.info.TCGA.RData",sep=""))

########################################################################

tumor.samples=sampleinfo[which(sampleinfo$sample_type=="Tumor"),]

########################################################################

## get tx2gene info

samples=tumor.samples$id

sinfo=read.table(paste(pathExpression, annot, "/",samples[1],"/abundance.tsv",sep=""), h=T, stringsAsFactors=F)
geneid=unlist(lapply(sinfo$target_id, function(x) unlist(strsplit(x, split=":"))[1]))

tx2gene=data.frame("txid"=sinfo$target_id, "geneid"=geneid)

########################################################################

## assemble kallisto counts

files=paste(pathExpression, "/", annot,"/", samples, "/abundance.h5", sep="")
names(files)=samples

txi.kallisto <- tximport(files, type = "kallisto", tx2gene = tx2gene)

########################################################################

## AJCC staging

ajcc=rep(NA, nrow(tumor.samples))

ajcc[which(tumor.samples$tumor_stage%in%c("stage i", "stage ii"))]="12"
ajcc[which(tumor.samples$tumor_stage%in%c("stage iii", "stage iiia", "stage iiic", "stage iv"))]="34"

colData=data.frame("Sex"=as.factor(tumor.samples$gender), "AJCC"=ajcc)

dds=DESeqDataSetFromTximport(txi.kallisto, colData=colData, design = ~Sex+AJCC)

dds=DESeq(dds, test="Wald",  minReplicatesForReplace=10, parallel=T)

########################################################################

results=results(dds, contrast=c("AJCC", "34", "12"))

results=lfcShrink(dds, coef="AJCC_34_vs_12", res=results, type="apeglm")

results=as.data.frame(results)

results=results[order(results$padj),]

########################################################################

write.table(results, file=paste(pathDifferentialExpression, annot, "/DifferentialExpression_AJCC_34_vs_12.txt",sep=""), row.names=T, col.names=T, sep="\t", quote=F)

########################################################################
