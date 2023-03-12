###########################################################################

pathDocs="../../docs/"
pathExpression="../../results/expression_estimation/"
pathSplicing="../../results/splicing_analysis/"

expdata="AllTranscripts_Ensembl109_noMT_norRNA_nohaplo"

###########################################################################

tumor.samples=read.table(paste(pathDocs, "TumorSamples_Ng2022.txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
tumor.samples=tumor.samples[which(tumor.samples$tumor_biopsyID!=""),]

###########################################################################

liver.samples=read.table(paste(pathDocs, "LiverSamples_Ng2022.txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
liver.samples=liver.samples[which(liver.samples$biopsyID!="POOL"),]

###########################################################################

sample.info.t=data.frame("BiopsyID"=tumor.samples$tumor_biopsyID, "TissueType"=rep("Tumor", nrow(tumor.samples)), "Sex"=tumor.samples$sex, "AgeAtBiopsy"=tumor.samples$age_at_biopsy, "EdmondsonGrade"=tumor.samples$edmondson, "Cirrhosis"=tumor.samples$cirrhosis, "Diseases"=tumor.samples$underlying_liver_disease)

sample.info.l=data.frame("BiopsyID"=liver.samples$biopsyID, "TissueType"=rep("Liver", nrow(liver.samples)), "Sex"=liver.samples$sex, "AgeAtBiopsy"=liver.samples$age_at_biopsy, "EdmondsonGrade"=rep(NA, nrow(liver.samples)), "Cirrhosis"=rep(NA, nrow(liver.samples)), "Diseases"=rep(NA, nrow(liver.samples)))

sample.info=rbind(sample.info.t, sample.info.l)

###########################################################################

tpm=read.table(paste(pathExpression, expdata, "/AllSamples_KallistoRawTPM.txt",sep=""),h=T, stringsAsFactors=F)


###########################################################################

## select one sample per patient for DE

dupli=which(duplicated(tumor.samples$Patient_ID))

if(length(dupli)>0){
    tumor.samples=tumor.samples[-dupli,]
}

###########################################################################


write.table(sample.info, paste(pathResults, "SampleInfo.txt", sep=""), row.names=F, col.names=T, sep="\t", quote=F)

###########################################################################
