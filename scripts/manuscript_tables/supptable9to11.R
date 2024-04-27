#############################################################################

if(!("pathTables"%in%objects())){
    pathRData="../../data_for_publication/RData/"
    pathTables="../../data_for_publication/tables/"

    library(vioplot)

    load=TRUE
    prepare=TRUE
}

#############################################################################

if(load){
    ## cited lncRNAs
    load(paste(pathRData, "data.PubMed.analysis.RData",sep=""))

    ## gene info
    load(paste(pathRData, "data.gene.info.RData", sep=""))

    ## differential expression
    load(paste(pathRData, "data.diffexp.RData",sep=""))

    load=FALSE
}

#############################################################################

if(prepare){
    for(type in c("tnt.meric", "tnt.tcga", "grades")){
        this.de=get(paste("diffexp",type,sep="."))
        this.de$GeneID=rownames(this.de)
        this.de$GeneName=geneinfo[rownames(this.de), "Name"]
        this.de$GeneBiotype=geneinfo[rownames(this.de), "Biotype"]

        this.de$NbCitations=rep(NA, nrow(this.de))
        this.de[pc, "NbCitations"]=rep(0, length(pc))
        this.de[lnc, "NbCitations"]=rep(0, length(lnc))
        this.de[names(nb.citations.pc), "NbCitations"]=nb.citations.pc
        this.de[names(nb.citations.lnc), "NbCitations"]=nb.citations.lnc

        cols=c("GeneID","GeneName", "GeneBiotype", "NbCitations")
        this.de=this.de[,c(cols, setdiff(colnames(this.de), cols))]

        assign(paste("diffexp",type,sep="."), this.de)
    }

    prepare=F
}

#############################################################################
## supp table 9 - tnt meric
write.table(diffexp.tnt.meric, paste(pathTables, "SupplementaryTable9.tsv",sep=""), sep="\t", quote=F, row.names=F, col.names=T)

## supp table 10 - grades meric
write.table(diffexp.grades, paste(pathTables, "SupplementaryTable10.tsv",sep=""), sep="\t", quote=F, row.names=F, col.names=T)

## supp table 11 - tnt tcga
write.table(diffexp.tnt.tcga, paste(pathTables, "SupplementaryTable11.tsv",sep=""), sep="\t", quote=F, row.names=F, col.names=T)

#############################################################################
