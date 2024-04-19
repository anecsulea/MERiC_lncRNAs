##########################################################################

if(!("pathFigures"%in%objects())){
    pathRData="../../data_for_publication/RData/"
    pathTables="../../data_for_publication/main_tables/"

    options(stringsAsFactors=F)

    load=TRUE
    prepare=TRUE
}

##########################################################################

if(load){
    ## info for each publication
    load(paste(pathRData, "data.PubMed.analysis.RData", sep=""))

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

    ## select those with at least 10 citations

    results.lnc=results.lnc[which(results.lnc$NbCitations>=10),]

    citations=unlist(lapply(results.lnc$PubMedIDs, function(x) unlist(strsplit(x, split=";"))))

    id=rep(results.lnc$ID, results.lnc$NbCitations)
    name=rep(results.lnc$Name, results.lnc$NbCitations)

    results=data.frame("ID"=id, "Name"=name, "PubMedID"=citations)

    prepare=FALSE
}

##########################################################################

## write raw table - we will add manual assessment later

write.table(results, file=paste(pathTables, "SupplementaryTable9_rawdata.tsv",sep=""), row.names=F, col.names=T, sep="\t", quote=F)

##########################################################################
