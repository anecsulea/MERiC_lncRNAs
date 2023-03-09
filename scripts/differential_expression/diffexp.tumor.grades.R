########################################################################

library(DESeq2)
library(BiocParallel)

register(MulticoreParam(4))

########################################################################

pathExpression="../../results/expression_estimation/"
pathDifferentialExpression="../../results/differential_expression/"
pathAnnot="../../data/ensembl_annotations/"

annot="AllTranscripts_Ensembl109_noMT_norRNA"

########################################################################

read.counts=read.table(paste(pathExpression, annot, "/AllSamples_KallistoEstimatedCounts.txt",sep=""), h=T, stringsAsFactors=F)

read.counts=round(read.counts)

########################################################################

sampleinfo=read.table(paste(pathDifferentialExpression, "SampleInfo.txt", sep=""), h=T, stringsAsFactors=F,sep="\t")

sampleinfo=sampleinfo[which(sampleinfo$BiopsyID%in%colnames(read.counts)),]

sampleinfo=sampleinfo[which(sampleinfo$TissueType=="Tumor"),]

read.counts=read.counts[,sampleinfo$BiopsyID]

########################################################################

geneinfo=read.table(paste(pathAnnot, "GeneInfo_Ensembl109.txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")

geneinfo=geneinfo[which(geneinfo$Gene.stable.ID%in%rownames(read.counts)),]

pc=geneinfo$Gene.stable.ID[which(geneinfo$Gene.type=="protein_coding")]

lnc=geneinfo$Gene.stable.ID[which(geneinfo$Gene.type=="lncRNA")]

read.counts=read.counts[c(pc, lnc),]

########################################################################

eg=rep(NA, nrow(sampleinfo))

eg[which(sampleinfo$EdmondsonGrade%in%c(1,2))]="12"
eg[which(sampleinfo$EdmondsonGrade%in%c(3,4))]="34"

colData=data.frame("Sex"=as.factor(sampleinfo$Sex), "EdmondsonGrade"=eg)

dds=DESeqDataSetFromMatrix(read.counts, colData=colData, design = ~Sex+EdmondsonGrade)

dds=DESeq(dds, test="Wald",  minReplicatesForReplace=10, parallel=T)

########################################################################

res.tissue.reduced=results(dds, contrast=c("EdmondsonGrade", "34", "12"))

res.tissue.reduced=lfcShrink(dds, coef="EdmondsonGrade_34_vs_12", res=res.tissue.reduced, type="apeglm")

res.tissue.reduced=as.data.frame(res.tissue.reduced)

res.tissue.reduced=res.tissue.reduced[order(res.tissue.reduced$padj),]

########################################################################

write.table(res.tissue.reduced, file=paste(pathDifferentialExpression, annot, "/DifferentialExpression_EdmondsonGrade_34_vs_12.txt",sep=""), row.names=T, col.names=T, sep="\t", quote=F)

########################################################################
