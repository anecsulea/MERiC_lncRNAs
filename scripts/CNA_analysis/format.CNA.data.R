######################################################################

pathData="../../data/CNA/2022/"
pathResults="../../results/CNA_analysis/"

######################################################################

cre=read.table(paste(pathData, "all_CRE.geneCN.ampdel.txt",sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
cre=cre[which(!is.na(cre$seg_start)),]

allex=read.table(paste(pathData, "all_AllExonv6.geneCN.ampdel.txt",sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
allex=allex[which(!is.na(allex$seg_start)),]

######################################################################

cre.bed=cre[,c("seg_chr", "seg_start", "seg_end")]
colnames(cre.bed)=c("chr", "start", "end")
cre.bed$id=paste(cre.bed$chr,":", cre.bed$start, "-", cre.bed$end,sep="")
cre.bed$chr=paste("chr", cre.bed$chr,sep="")

write.table(cre.bed, file=paste(pathResults, "CRE_ampdel_hg19.bed",sep=""), row.names=F, col.names=F, sep="\t", quote=F)

######################################################################

allex.bed=allex[,c("seg_chr", "seg_start", "seg_end")]
colnames(allex.bed)=c("chr", "start", "end")
allex.bed$id=paste(allex.bed$chr,":", allex.bed$start, "-", allex.bed$end,sep="")
allex.bed$chr=paste("chr", allex.bed$chr,sep="")

write.table(allex.bed, file=paste(pathResults, "AllExonsv6_ampdel_hg19.bed",sep=""), row.names=F, col.names=F, sep="\t", quote=F)

######################################################################
