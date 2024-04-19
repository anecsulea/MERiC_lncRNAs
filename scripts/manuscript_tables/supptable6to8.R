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

    ## gene info
    load(paste(pathRData, "data.gene.info.RData", sep=""))

    ## TPM MERiC
    load(paste(pathRData, "data.expression.levels.MERiC.RData",sep=""))
    tpm.meric=tpm

    ## sample info
    load(paste(pathRData, "data.sample.info.MERiC.RData",sep=""))
    tumor.samples.meric=tumor.samples
    nontumor.samples.meric=nontumor.samples

    load=FALSE
}

##########################################################################

if(prepare){

    tpm.meric=tpm.meric[c(pc, lnc), c(tumor.samples.meric$tumor_biopsyID, nontumor.samples.meric$biopsyID)]

    prepare=FALSE
}

##########################################################################

## info for selected tumor samples

write.table(tumor.samples.meric,  file=paste(pathTables, "SupplementaryTable5.tsv",sep=""), row.names=F, col.names=T, sep="\t", quote=F)

## info for selected adjacent tissue samples

write.table(nontumor.samples.meric,  file=paste(pathTables, "SupplementaryTable6.tsv",sep=""), row.names=F, col.names=T, sep="\t", quote=F)

## TPM values for pc and lncRNA genes

write.table(tpm.meric[,tumor.samples.meric$tumor_biopsyID], file=paste(pathTables, "SupplementaryTable8.tsv",sep=""), row.names=T, col.names=T, sep="\t", quote=F)


##########################################################################
