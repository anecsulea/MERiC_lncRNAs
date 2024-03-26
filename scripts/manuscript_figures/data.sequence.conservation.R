##########################################################################

pathConservation="../../results/sequence_conservation/"
pathRData="../../data_for_publication/RData/"

##########################################################################

for(phast in c("30way", "100way")){
    
    phastcons=read.table(paste(pathConservation, phast,"/PhastCons_GeneAverage_EnsemblNoOverlaps.txt", sep=""), h=T, stringsAsFactors=F)

    rownames(phastcons)=phastcons$Gene

    save(phastcons, file=paste(pathRData, "data.phastcons.",phast,".RData",sep=""))
}

##########################################################################
