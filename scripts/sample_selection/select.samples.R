###########################################################################

pathDocs="../../docs/"
pathExpression="../../results/expression_estimation/"
pathSplicing="../../results/splicing_analysis/"

expdata="AllTranscripts_Ensembl109_noMT_norRNA_nohaplo"
annot="AllTranscripts_Ensembl109"

###########################################################################

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

tpm=read.table(paste(pathExpression, expdata, "/AllSamples_KallistoRawTPM.txt",sep=""),h=T, stringsAsFactors=F)

nb.exp.genes=apply(tpm,2,function(x) length(which(x>=1)))
names(nb.exp.genes)=colnames(tpm)

###########################################################################

ir=read.table(paste(pathSplicing, "/AllSamples_IntronExonRatio_", annot, ".txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
rownames(ir)=ir$GeneID
ir=ir[,which(colnames(ir)!="GeneID")]

median.ir=apply(ir,2, median, na.rm=T)
names(median.ir)=colnames(ir)

###########################################################################

sample.info$NbExpressedGenes=nb.exp.genes[sample.info$BiopsyID]
sample.info$MedianIntronExonRatio=median.ir[sample.info$BiopsyID]

###########################################################################

sample.info$Selected=rep("Yes", nrow(sample.info))
sample.info$Selected[which(sample.info$MedianIntronExonRatio>=0.15)]="No"
sample.info$Selected[which(is.na(sample.info$NbExpressedGenes))]="No"

###########################################################################

## select one sample per patient for DE - highest grade if several

sel=sample.info[which(sample.info$Selected=="Yes"),]
sel=sel[order(sel$NbExpressedGenes, decreasing=T),]
sel=sel[order(sel$EdmondsonGrade, decreasing=T),]
dupli=which(duplicated(sel$PatientID))

if(length(dupli)>0){
    print("keeping one sample per patient")
    id.dupli=sel$BiopsyID[dupli]
    sample.info[which(sample.info$BiopsyID%in%id.dupli), "Selected"]="No"
}

###########################################################################

write.table(sample.info, paste(pathResults, "AllSampleInfo.txt", sep=""), row.names=F, col.names=T, sep="\t", quote=F)

###########################################################################

selected.info=sample.info[which(sample.info$Selected=="Yes"),]
selected.info=selected.info[,which(!colnames(selected.info)%in%c("Selected"))]
selected.info$PatientID[which(selected.info$TissueType=="Liver")]=NA

write.table(selected.info, paste(pathResults, "AnalyzedSamples.txt", sep=""), row.names=F, col.names=T, sep="\t", quote=F)

###########################################################################
