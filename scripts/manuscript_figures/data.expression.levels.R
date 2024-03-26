##########################################################################

pathExpression="../../results/expression_estimation/"
pathExpressionTCGA="../../results/expression_estimation_TCGA/"
pathRData="../../data_for_publication/RData/"

annot="AllTranscripts_Ensembl109_noMT_norRNA_nohaplo"

##########################################################################

tpm=read.table(paste(pathExpression, annot, "/AllSamples_KallistoRawTPM.txt",sep=""), h=T, stringsAsFactors=F, sep="\t")

##########################################################################

save(tpm, file=paste(pathRData, "data.expression.levels.MERiC.RData",sep=""))

##########################################################################

tpm=read.table(paste(pathExpressionTCGA, annot, "/AllSamples_KallistoRawTPM.txt",sep=""), h=T, stringsAsFactors=F, sep="\t")

##########################################################################

save(tpm, file=paste(pathRData, "data.expression.levels.TCGA.RData",sep=""))

##########################################################################
