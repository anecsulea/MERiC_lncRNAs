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

    load=FALSE
}

##########################################################################

if(prepare){
    ## for each lncRNA: how many citations, and which citations exactly

    unlisted.cited.lnc=unlist(lapply(pubs$CitedLnc, function(x) setdiff(unlist(strsplit(x, split=";")), "")))
    nb.cited.lnc=unlist(lapply(pubs$CitedLnc, function(x) length(setdiff(unlist(strsplit(x, split=";")), ""))))
    pmids.cited.lnc=rep(pubs$PMID, nb.cited.lnc)

    cited.pmids.per.lnc=tapply(pmids.cited.lnc, as.factor(unlisted.cited.lnc), function(x) paste(x, collapse=";"))

    nb.citations.per.lnc=as.numeric(table(unlisted.cited.lnc))
    names(nb.citations.per.lnc)=levels(as.factor(unlisted.cited.lnc))

    results.lnc=data.frame("ID"=intersect(unique(unlisted.cited.lnc), lnc))
    results.lnc$Name=geneinfo[results.lnc$ID, "Name"]
    results.lnc$NbCitations=nb.citations.per.lnc[results.lnc$ID]
    results.lnc$PubMedIDs=cited.pmids.per.lnc[results.lnc$ID]

    results.lnc=results.lnc[order(results.lnc$NbCitations, decreasing=T),]


    ## for each protein-coding gene: how many citations, and which citations exactly

    unlisted.cited.pc=unlist(lapply(pubs$CitedPc, function(x) setdiff(unlist(strsplit(x, split=";")), "")))
    nb.cited.pc=unlist(lapply(pubs$CitedPc, function(x) length(setdiff(unlist(strsplit(x, split=";")), ""))))
    pmids.cited.pc=rep(pubs$PMID, nb.cited.pc)

    cited.pmids.per.pc=tapply(pmids.cited.pc, as.factor(unlisted.cited.pc), function(x) paste(x, collapse=";"))

    nb.citations.per.pc=as.numeric(table(unlisted.cited.pc))
    names(nb.citations.per.pc)=levels(as.factor(unlisted.cited.pc))

    results.pc=data.frame("ID"=intersect(unique(unlisted.cited.pc), pc))
    results.pc$Name=geneinfo[results.pc$ID, "Name"]
    results.pc$NbCitations=nb.citations.per.pc[results.pc$ID]
    results.pc$PubMedIDs=cited.pmids.per.pc[results.pc$ID]

    results.pc=results.pc[order(results.pc$NbCitations, decreasing=T),]

    prepare=FALSE
}

##########################################################################

## list of all publications

write.table(pubs, file=paste(pathTables, "SupplementaryTable1.tsv",sep=""), row.names=F, col.names=T, sep="\t", quote=F)

##########################################################################

## list of cited lncRNAs

write.table(results.lnc, file=paste(pathTables, "SupplementaryTable2.tsv",sep=""), row.names=F, col.names=T, sep="\t", quote=F)


## list of cited protein-coding genes

write.table(results.pc, file=paste(pathTables, "SupplementaryTable3.tsv",sep=""), row.names=F, col.names=T, sep="\t", quote=F)

##########################################################################
