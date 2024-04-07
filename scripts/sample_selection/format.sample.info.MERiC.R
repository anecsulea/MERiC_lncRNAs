###########################################################################

pathDocs="../../docs/"
pathRData="../../data_for_publication/RData/"

###########################################################################

sample.annot=read.table(paste(pathDocs, "SampleAnnotation_CN_AN.txt",sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
rownames(sample.annot)=sample.annot$biopsy_id
sample.annot$Patient_ID=as.character(sample.annot$Patient_ID)

sample.annot=sample.annot[which(sample.annot$comment.AN==""),]

###########################################################################

tumor.samples=read.table(paste(pathDocs, "TumorSamples_Ng2022.txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
tumor.samples=tumor.samples[which(tumor.samples$tumor_biopsyID!=""),]

selected.columns=c("tumor_biopsyID" , "Patient_ID", "sex", "age_at_biopsy", "bclc", "cirrhosis", "child", "meld", "underlying_liver_disease", "number.of.tumors", "macro_vascular_invasion", "metastasis", "tumor_location", "edmondson", "percent_tumor", "immunophenotype")

tumor.samples=tumor.samples[,selected.columns]

tumor.samples$Patient_ID=as.character(tumor.samples$Patient_ID)

patient.sex=tumor.samples$sex
names(patient.sex)=tumor.samples$Patient_ID

###########################################################################

tumor.samples=tumor.samples[which(tumor.samples$tumor_biopsyID%in%sample.annot$biopsy_id),]

patients=sample.annot[tumor.samples$tumor_biopsyID, "Patient_ID"]

annot.nontumor.samples=sample.annot[which(sample.annot$Patient_ID%in%patients & sample.annot$tissue=="Non-Tumor-Liver"),]

###########################################################################

tumor.samples.unique=tumor.samples[which(!duplicated(tumor.samples$Patient_ID)),]

patient.synonyms=tumor.samples.unique$Patient_ID

names(patient.synonyms)=sample.annot[tumor.samples.unique$tumor_biopsyID, "Patient_ID"]

nontumor.samples=annot.nontumor.samples[,c("biopsy_id", "Patient_ID", "bio_age", "diagnosis1", "diagnosis2", "diagnosis3")]
colnames(nontumor.samples)=c("biopsyID", "Patient_ID", "age_at_biopsy",  "diagnosis1", "diagnosis2", "diagnosis3")
nontumor.samples$Patient_ID=patient.synonyms[nontumor.samples$Patient_ID]

nontumor.samples$sex=patient.sex[nontumor.samples$Patient_ID]

###########################################################################

save(list=c("tumor.samples","nontumor.samples"), file=paste(pathRData, "data.sample.info.MERiC.RData",sep=""))

###########################################################################

write.table(tumor.samples, paste(pathDocs, "MERiC_Samples_Tumor.txt",sep=""), row.names=F, col.names=T, sep="\t", quote="\"")

write.table(nontumor.samples, paste(pathDocs, "MERiC_Samples_AdjacentTissue.txt",sep=""), row.names=F, col.names=T, sep="\t", quote="\"")

###########################################################################

