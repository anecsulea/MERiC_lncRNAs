##########################################################################

if(!("pathFigures"%in%objects())){
    pathRData="../../data_for_publication/RData/"
    pathTables="../../data_for_publication/tables/"

    options(stringsAsFactors=F)

    load=TRUE
    prepare=TRUE
}

##########################################################################

if(load){
    ## info for each publication
    load(paste(pathRData, "data.PubMed.analysis.2023.RData", sep=""))

    ## gene info
    load(paste(pathRData, "data.gene.info.RData", sep=""))

    ## gene overlaps
    load(paste(pathRData, "data.gene.overlaps.RData", sep=""))

    ## TPM MERiC
    load(paste(pathRData, "data.expression.levels.MERiC.RData",sep=""))
    tpm.meric=tpm

    ## TPM TCGA
    load(paste(pathRData, "data.expression.levels.TCGA.RData",sep=""))
    tpm.tcga=tpm

    ## sample info
    load(paste(pathRData, "data.sample.info.MERiC.RData",sep=""))
    tumor.samples.meric=tumor.samples
    nontumor.samples.meric=nontumor.samples

    ## sample info TCGA

    load(paste(pathRData, "data.sample.info.TCGA.RData",sep=""))
    sampleinfo.tcga=sampleinfo
    tumor.samples.tcga=sampleinfo.tcga[which(sampleinfo$sample_type=="Tumor"),]
    nontumor.samples.tcga=sampleinfo.tcga[which(sampleinfo$sample_type=="Non-Tumor"),]

    ## sequence conservation
    load(paste(pathRData, "data.phastcons.30way.RData",sep=""))

    load=FALSE
}

##########################################################################

if(prepare){

    ## extract overlaps with pc genes
    overlaps.pc=gene.overlaps[which(gene.overlaps$NeighborID%in%pc),]

    sense.overlaps.pc=overlaps.pc[which(overlaps.pc$Type=="sense"),]
    antisense.overlaps.pc=overlaps.pc[which(overlaps.pc$Type=="antisense"),]

    sense.overlaps=gene.overlaps[which(gene.overlaps$Type=="sense"),]
    antisense.overlaps=gene.overlaps[which(gene.overlaps$Type=="antisense"),]

    sense.overlapping.genes=tapply(sense.overlaps$NeighborID, as.factor(sense.overlaps$GeneID), function(x) paste(unique(x), collapse=";"))
    antisense.overlapping.genes=tapply(antisense.overlaps$NeighborID, as.factor(antisense.overlaps$GeneID), function(x) paste(unique(x), collapse=";"))

    ## bidirectional promoters with pc genes

    biprom1kb=biprom1kb[which(!is.na(biprom1kb$GenesCloseTSS)),]
    biprom5kb=biprom5kb[which(!is.na(biprom5kb$GenesCloseTSS)),]

    biprom1kb.genes=tapply(biprom1kb$GenesCloseTSS, as.factor(biprom1kb$GeneID), function(x) paste(unique(unlist(lapply(x, function(y) unlist(strsplit(y, split=","))))), collapse=","))
    biprom5kb.genes=tapply(biprom5kb$GenesCloseTSS, as.factor(biprom5kb$GeneID), function(x) paste(unique(unlist(lapply(x, function(y) unlist(strsplit(y, split=","))))), collapse=","))

    ## average expression level across various samples
    meantpm.nontumor.meric=apply(tpm.meric[ ,nontumor.samples.meric$biopsyID],1, mean)
    meantpm.tumor.meric=apply(tpm.meric[ ,tumor.samples.meric$tumor_biopsyID],1, mean)

    meantpm.tumor.tcga=apply(tpm.tcga[, tumor.samples.tcga$id], 1, mean)
    meantpm.nontumor.tcga=apply(tpm.tcga[, nontumor.samples.tcga$id], 1, mean)

    ## prepare table with info for all lncRNAs

    info.lnc=geneinfo[lnc,]
    rownames(info.lnc)=info.lnc$ID

    info.lnc$NbCitations=rep(0, nrow(info.lnc))
    info.lnc[names(nb.citations.lnc), "NbCitations"]=nb.citations.lnc
    info.lnc$SenseOverlaps=sense.overlapping.genes[info.lnc$ID]
    info.lnc$AntisenseOverlaps=antisense.overlapping.genes[info.lnc$ID]

    info.lnc$BidirectionalPromoters1kb=biprom1kb.genes[info.lnc$ID]
    info.lnc$BidirectionalPromoters5kb=biprom5kb.genes[info.lnc$ID]

    info.lnc$PhastConsScore=phastcons[info.lnc$ID, "Score"]

    info.lnc$MeanTPM.AdjacentTissue.MERiC=meantpm.nontumor.meric[info.lnc$ID]
    info.lnc$MeanTPM.AdjacentTissue.TCGA=meantpm.nontumor.tcga[info.lnc$ID]

    info.lnc$MeanTPM.Tumor.MERiC=meantpm.tumor.meric[info.lnc$ID]
    info.lnc$MeanTPM.Tumor.TCGA=meantpm.tumor.tcga[info.lnc$ID]

     ## prepare table with info for all protein-coding genes

    info.pc=geneinfo[pc,]
    rownames(info.pc)=info.pc$ID

    info.pc$NbCitations=rep(0, nrow(info.pc))
    info.pc[names(nb.citations.pc), "NbCitations"]=nb.citations.pc
    info.pc$SenseOverlaps=sense.overlapping.genes[info.pc$ID]
    info.pc$AntisenseOverlaps=antisense.overlapping.genes[info.pc$ID]

    info.pc$BidirectionalPromoters1kb=biprom1kb.genes[info.pc$ID]
    info.pc$BidirectionalPromoters5kb=biprom5kb.genes[info.pc$ID]

    info.pc$PhastConsScore=phastcons[info.pc$ID, "Score"]

    info.pc$MeanTPM.AdjacentTissue.MERiC=meantpm.nontumor.meric[info.pc$ID]
    info.pc$MeanTPM.AdjacentTissue.TCGA=meantpm.nontumor.tcga[info.pc$ID]

    info.pc$MeanTPM.Tumor.MERiC=meantpm.tumor.meric[info.pc$ID]
    info.pc$MeanTPM.Tumor.TCGA=meantpm.tumor.tcga[info.pc$ID]

}

##########################################################################

## info for lncRNAs

write.table(info.lnc, file=paste(pathTables, "SupplementaryTable6.tsv",sep=""), row.names=F, col.names=T, sep="\t", quote=F)

## info for protein-coding genes

write.table(info.pc, file=paste(pathTables, "SupplementaryTable7.tsv",sep=""), row.names=F, col.names=T, sep="\t", quote=F)

##########################################################################
