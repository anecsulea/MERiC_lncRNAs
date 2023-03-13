########################################################################

library(DESeq2)
library(BiocParallel)

register(MulticoreParam(4))

########################################################################

pathExpression="../../results/expression_estimation/"
pathDifferentialExpression="../../results/differential_expression/"
pathAnnot="../../data/ensembl_annotations/"
pathSampleInfo="../../results/sample_info/"

annot="AllTranscripts_Ensembl109_noMT_norRNA_nohaplo"

########################################################################

read.counts=read.table(paste(pathExpression, annot, "/AllSamples_KallistoEstimatedCounts.txt",sep=""), h=T, stringsAsFactors=F)

read.counts=round(read.counts)

########################################################################

sampleinfo=read.table(paste(pathSampleInfo, "AnalyzedSamples.txt", sep=""), h=T, stringsAsFactors=F,sep="\t")

print(paste("have reads for all samples: ",all(sampleinfo$BiopsyID%in%colnames(read.counts))))

sampleinfo=sampleinfo[which(sampleinfo$TissueType=="Tumor"),]

read.counts=read.counts[,sampleinfo$BiopsyID]

########################################################################

geneinfo=read.table(paste(pathAnnot, "GeneInfo_Ensembl109.txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")

geneinfo=geneinfo[which(geneinfo$Gene.stable.ID%in%rownames(read.counts)),]

pc=geneinfo$Gene.stable.ID[which(geneinfo$Gene.type=="protein_coding")]

lnc=geneinfo$Gene.stable.ID[which(geneinfo$Gene.type=="lncRNA")]

pseudo=geneinfo$Gene.stable.ID[which(geneinfo$Gene.type%in%c("transcribed_unitary_pseudogene", "transcribed_unprocessed_pseudogene"))]

read.counts=read.counts[c(pc, lnc, pseudo),]

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
