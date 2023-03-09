########################################################################

pathExpression="../../results/expression_estimation/"
pathDifferentialExpression="../../results/differential_expression/"
pathAnnot="../../data/ensembl_annotations/"

annot="AllTranscripts_Ensembl109_noMT_norRNA"

########################################################################

library(ade4)

########################################################################

tpm=read.table(paste(pathExpression, annot, "/AllSamples_KallistoRawTPM.txt",sep=""), h=T, stringsAsFactors=F)

########################################################################

sampleinfo=read.table(paste(pathDifferentialExpression, "SampleInfo.txt", sep=""), h=T, stringsAsFactors=F,sep="\t")

col.samples=rep("gray40", nrow(sampleinfo))
col.samples[which(sampleinfo$TissueType=="Tumor" & sampleinfo$EdmondsonGrade==1)]="yellow"
col.samples[which(sampleinfo$TissueType=="Tumor" & sampleinfo$EdmondsonGrade==2)]="orange"
col.samples[which(sampleinfo$TissueType=="Tumor" & sampleinfo$EdmondsonGrade==3)]="darkorange"
col.samples[which(sampleinfo$TissueType=="Tumor" & sampleinfo$EdmondsonGrade==4)]="red"

names(col.samples)=sampleinfo$BiopsyID

########################################################################

geneinfo=read.table(paste(pathAnnot, "GeneInfo_Ensembl109.txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")

pc=geneinfo$Gene.stable.ID[which(geneinfo$Gene.type=="protein_coding")]

lnc=geneinfo$Gene.stable.ID[which(geneinfo$Gene.type=="lncRNA")]

########################################################################

tpm.pc=tpm[intersect(pc, rownames(tpm)),]

tpm.lnc=tpm[intersect(lnc, rownames(tpm)),]

########################################################################

pca.pc=dudi.pca(t(log2(tpm.pc+1)), scale=F, center=T, scannf=F, nf=5)

pca.lnc=dudi.pca(t(log2(tpm.lnc+1)), scale=F, center=T, scannf=F, nf=5)

########################################################################

par(mfrow=c(1,2))

plot(pca.pc$li[,1], pca.pc$li[,2], pch=20, col=col.samples[rownames(pca.pc$li)])

plot(pca.lnc$li[,1], pca.lnc$li[,2], pch=20, col=col.samples[rownames(pca.lnc$li)])

########################################################################
