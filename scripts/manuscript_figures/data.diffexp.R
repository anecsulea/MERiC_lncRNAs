##########################################################################

pathDifferentialExpression="../../results/differential_expression/"
pathDifferentialExpressionTCGA="../../results/differential_expression_TCGA/"
pathRData="../../data_for_publication/RData/"

annot="AllTranscripts_Ensembl109_noMT_norRNA_nohaplo"

##########################################################################

## MERIC

diffexp.grades=read.table(paste(pathDifferentialExpression, annot, "/DifferentialExpression_EdmondsonGrade_34_vs_12.txt", sep=""), h=T, stringsAsFactors=F, sep="\t")

diffexp.tissues=read.table(paste(pathDifferentialExpression, annot, "/DifferentialExpression_Tumor_vs_Liver.txt", sep=""), h=T, stringsAsFactors=F, sep="\t")

## TCGA

diffexp.tnt=read.table(paste(pathDifferentialExpressionTCGA, annot, "/DifferentialExpression_Tumor_vs_NonTumor.txt", sep=""), h=T, stringsAsFactors=F, sep="\t")

##########################################################################

save(list=c("diffexp.grades", "diffexp.tissues", "diffexp.tnt"), file=paste(pathRData, "data.diffexp.RData",sep=""))

##########################################################################
