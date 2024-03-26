##########################################################################

if(!("pathFigures"%in%objects())){
    pathRData="../../data_for_publication/RData/"
    pathFigures="../../data_for_publication/main_figures/"

    load=TRUE
    prepare=TRUE
}

##########################################################################

if(load){
    ## cited lncRNAs
    load(paste(pathRData, "data.PubMed.analysis.RData", sep=""))

    ## gene info
    load(paste(pathRData, "data.gene.info.RData", sep=""))

    ## gene overlaps
    load(paste(pathRData, "data.gene.overlaps.RData", sep=""))

    ## TPM
    load(paste(pathRData, "data.expression.levels.MERiC.RData",sep=""))
    tpm.meric=tpm

    load(paste(pathRData, "data.expression.levels.TCGA.RData",sep=""))
    tpm.tcga=tpm    

    load=FALSE
}

##########################################################################

if(prepare){
    ## other lnc, not cited
    other.lnc=setdiff(lnc, all.cited.lnc)

    ## extract overlaps with pc genes
    overlaps.pc=gene.overlaps[which(gene.overlaps$NeighborID%in%pc),]

    ## sense and antisense overlaps
    sense.overlaps.pc=overlaps.pc[which(overlaps.pc$Type=="sense"),]
    antisense.overlaps.pc=overlaps.pc[which(overlaps.pc$Type=="antisense"),]
    
    ## bidirectional promoters with pc genes

    biprom1pc=biprom1kb[which(biprom1kb$GenesCloseTSS%in%pc),]
    biprom5pc=biprom1kb[which(biprom5kb$GenesCloseTSS%in%pc),]

    ## 
    
    prepare=FALSE
}

##########################################################################

pdf(paste(pathFigures, "Figure2.pdf", sep=""), width=6.85, height=3.5)

##########################################################################


##########################################################################

##########################################################################

dev.off()

##########################################################################
