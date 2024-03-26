##########################################################################

pathConservation="../../results/sequence_conservation/"
pathRData="../../data_for_publication/RData/"

##########################################################################

phastcons=read.table(paste(pathConservation, "30way/PhastCons_GeneAverage_EnsemblNoOverlaps.txt", sep=""), h=T, stringsAsFactors=F)

rownames(phastcons)=phastcons$Gene

##########################################################################

save(phastcons, file=paste(pathRData, "data.sequence.conservation.RData",sep=""))

##########################################################################
