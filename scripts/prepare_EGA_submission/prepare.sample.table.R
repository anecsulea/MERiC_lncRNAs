#############################################################################

pathTables="../../data_for_publication/tables/"
pathEGA="../../data_for_publication/EGA_submission/"

#############################################################################

sampleinfo=read.table(paste(pathTables, "SupplementaryTable5.tsv", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
rownames(sampleinfo)=sampleinfo$biopsyID

#############################################################################

columnOrder=c("title","alias","description","subjectId","bioSampleId","caseOrControl","gender","organismPart","cellLine","region","phenotype")

#############################################################################

results=data.frame("title"=sampleinfo$biopsyID, "alias"=sampleinfo$biopsyID, "description"=rep("RNA-seq from adjacent tissue biopsy", dim(sampleinfo)[1]),  "subjectId"=sampleinfo$Patient_ID, "bioSampleId"=rep("", dim(sampleinfo)[1]), "caseOrControl"=rep("case", dim(sampleinfo)[1]), "gender"=sampleinfo$sex, "organismPart"=rep("liver", dim(sampleinfo)[1]), "cellLine"=rep("", dim(sampleinfo)[1]), "region"=rep("adjacent tissue", dim(sampleinfo)[1]), "phenotype"=rep("hepatocellular carcinoma", dim(sampleinfo)[1]),  stringsAsFactors=F)

results$alias=paste("MERIC_",results$alias, sep="")

results$region[which(sampleinfo[results$title, "Tissue"]=="Tumor")]="tumor"
results$region[which(sampleinfo[results$title, "Tissue"]=="Non-Tumor")]="adjacent tissue"

results$gender[which(sampleinfo[results$title, "sex"]=="m")]="male"
results$gender[which(sampleinfo[results$title, "sex"]=="f")]="female"

#############################################################################

results=results[, columnOrder]

#############################################################################

write.table(results, file=paste(pathEGA, "EGA_submission_sampleinfo.csv", sep=""), row.names=F, col.names=T, sep=",", quote=T)

#############################################################################
