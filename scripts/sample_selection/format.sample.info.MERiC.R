###########################################################################

pathDocs="../../docs/"
pathRData="../../data_for_publication/RData/"

###########################################################################

sample.annot=read.table(paste(pathDocs, "SampleAnnotation_CN_AN.txt",sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")

###########################################################################

tumor.samples=read.table(paste(pathDocs, "TumorSamples_Ng2022.txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
tumor.samples=tumor.samples[which(tumor.samples$tumor_biopsyID!=""),]

selected.columns=c("tumor_biopsyID" , "Patient_ID", "sex", "age_at_biopsy", "bclc", "cirrhosis", "child", "meld", "underlying_liver_disease", "number.of.tumors", "macro_vascular_invasion", "metastasis", "tumor_location", "edmondson", "percent_tumor", "immunophenotype")

tumor.samples=tumor.samples[,selected.columns]

###########################################################################

liver.samples=read.table(paste(pathDocs, "LiverSamples_Ng2022.txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
liver.samples=liver.samples[which(liver.samples$biopsyID!="POOL"),]
liver.samples=liver.samples[which(liver.samples$RNAseq=="Y"),]

selected.columns=c("biopsyID", "sex", "age_at_biopsy", "normal_liver_quality")

liver.samples=liver.samples[,selected.columns]

###########################################################################

save(list=c("tumor.samples", "liver.samples"), file=paste(pathRData, "data.sample.info.RData",sep=""))

###########################################################################

