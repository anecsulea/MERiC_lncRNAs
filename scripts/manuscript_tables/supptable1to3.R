##########################################################################

if(!("pathFigures"%in%objects())){
    pathRData="../../data_for_publication/RData/"
    pathTables="../../data_for_publication/main_tables/"

    load=TRUE
    prepare=TRUE
}

##########################################################################

if(load){
    load(paste(pathRData, "data.PubMed.analysis.RData", sep=""))

    load=FALSE
}

##########################################################################

write.table(pubs, file=paste(pathTables, "SupplementaryTable1.tsv",sep=""), row.names=F, col.names=T, sep="\t", quote=F)

##########################################################################
