##########################################################################

pathDifferentialExpression="../../results/differential_expression/"
pathRData="../../data_for_publication/RData/"

annot="AllTranscripts_Ensembl109_noMT_norRNA_nohaplo"

##########################################################################

diffexp.grades=read.table(paste(pathDifferentialExpression, "DifferentialExpression_EdmondsonGrade_34_vs_12.txt", sep=""), h=T, stringsAsFactors=F, sep="\t")

diffexp.tissues=read.table(paste(pathDifferentialExpression, "DifferentialExpression_Tumor_vs_Liver.txt", sep=""), h=T, stringsAsFactors=F, sep="\t")

##########################################################################

save(list=c("diffexp.grades", "diffexp.tissues"), file=paste(pathRData, "data.diffexp.RData",sep=""))

##########################################################################
