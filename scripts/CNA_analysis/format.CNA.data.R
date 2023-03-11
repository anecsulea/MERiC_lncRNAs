######################################################################

pathData="../../data/CNA/2022/"
pathResults="../../results/CNA_analysis/"

######################################################################

cre=read.table(paste(pathData, "../../data/CNA/2022/all_CRE.geneCN.ampdel.txt",sep=""), h=T, stringsAsFactors=F)
allex=read.table(paste(pathData, "../../data/CNA/2022/all_AllExonv6.geneCN.ampdel.txt",sep=""), h=T, stringsAsFactors=F)

######################################################################
