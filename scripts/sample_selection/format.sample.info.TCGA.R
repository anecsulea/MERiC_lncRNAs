###########################################################################

pathTCGA="../../data/TCGA/"
pathRData="../../data_for_publication/RData/"

###########################################################################

## read sample info

sampleinfo=read.table(paste(pathTCGA, "all_info_paired_samples.txt", sep=""), h=T, stringsAsFactors=F, sep="\t")
sampleinfo$case_id=as.character(sampleinfo$case_id)
rownames(sampleinfo)=sampleinfo$id

sampleinfo$sample_type[which(sampleinfo$sample_type=="Solid Tissue Normal")]="Non-Tumor"
sampleinfo$sample_type[which(sampleinfo$sample_type=="Primary Tumor")]="Tumor"

########################################################################

save(list=c("sampleinfo"), file=paste(pathRData, "data.sample.info.TCGA.RData",sep=""))

###########################################################################

