##########################################################################

pathEnsembl="../../data/ensembl_annotations/"
pathRData="../../data_for_publication/RData/"

release=103

##########################################################################

geneinfo=read.table(paste(pathEnsembl, "GeneInfo_Ensembl",release, ".txt", sep=""), h=T, sep="\t", stringsAsFactors=F, quote="\"")

colnames(geneinfo)=c("ID", "Name", "Chr", "Start", "End", "Strand", "Biotype")
rownames(geneinfo)=geneinfo$ID

###############################################################################

lnc=geneinfo$ID[which(geneinfo$Biotype=="lncRNA")]
pc=geneinfo$ID[which(geneinfo$Biotype=="protein_coding")]

###############################################################################

save(list=c("geneinfo", "lnc", "pc"), file=paste(pathRData, "data.gene.info.RData",sep=""))

###############################################################################
