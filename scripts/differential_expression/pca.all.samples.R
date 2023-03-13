########################################################################

pathExpression="../../results/expression_estimation/"
pathDifferentialExpression="../../results/differential_expression/"
pathAnnot="../../data/ensembl_annotations/"
pathSampleInfo="../../results/sample_info/"
pathResults="../../results/figures/analyses/"

annot="AllTranscripts_Ensembl109_noMT_norRNA_nohaplo"

########################################################################

library(ade4)

########################################################################

tpm=read.table(paste(pathExpression, annot, "/AllSamples_KallistoRawTPM.txt",sep=""), h=T, stringsAsFactors=F)

########################################################################

tumor.samples=read.table(paste(pathDocs, "TumorSamples_Ng2022.txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
tumor.samples=tumor.samples[which(tumor.samples$tumor_biopsyID!=""),]

###########################################################################

liver.samples=read.table(paste(pathDocs, "LiverSamples_Ng2022.txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
liver.samples=liver.samples[which(liver.samples$biopsyID!="POOL"),]
liver.samples$PatientID=liver.samples$biopsyID

###########################################################################

sample.info.t=data.frame("BiopsyID"=tumor.samples$tumor_biopsyID, "TissueType"=rep("Tumor", nrow(tumor.samples)), "Sex"=tumor.samples$sex, "AgeAtBiopsy"=tumor.samples$age_at_biopsy, "EdmondsonGrade"=tumor.samples$edmondson, "Cirrhosis"=tumor.samples$cirrhosis, "Diseases"=tumor.samples$underlying_liver_disease, "PatientID"=tumor.samples$Patient_ID)

sample.info.l=data.frame("BiopsyID"=liver.samples$biopsyID, "TissueType"=rep("Liver", nrow(liver.samples)), "Sex"=liver.samples$sex, "AgeAtBiopsy"=liver.samples$age_at_biopsy, "EdmondsonGrade"=rep(NA, nrow(liver.samples)), "Cirrhosis"=rep(NA, nrow(liver.samples)), "Diseases"=rep(NA, nrow(liver.samples)), "PatientID"=liver.samples$PatientID)

sample.info=rbind(sample.info.t, sample.info.l)

###########################################################################

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

pdf(file=paste(pathResults, "PCA_all_samples_pc_lnc.pdf",sep=""),width=8,height=5)

par(mfrow=c(1,2))
par(mar=c(4.1,2.1,1.1,1.1))

plot(pca.pc$li[,1], pca.pc$li[,2], pch=20, col=col.samples[rownames(pca.pc$li)])

plot(pca.lnc$li[,1], pca.lnc$li[,2], pch=20, col=col.samples[rownames(pca.lnc$li)])

dev.off()

########################################################################
