##########################################################################

if(!("pathFigures"%in%objects())){
    pathRData="../../data_for_publication/RData/"
    pathFigures="../../data_for_publication/main_figures/"

    load=TRUE
    prepare=TRUE
}

##########################################################################

if(load){
    load(paste(pathRData, "data.PubMed.analysis.RData", sep=""))

    load=FALSE
}

##########################################################################
