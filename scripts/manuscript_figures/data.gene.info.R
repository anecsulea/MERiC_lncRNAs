##########################################################################

pathEnsembl="../../data/ensembl_annotations/"
pathRData="../../data_for_publication/RData/"

release=109

##########################################################################

geneinfo=read.table(paste(pathEnsembl, "GeneInfo_Ensembl",release, ".txt", sep=""), h=T, sep="\t", stringsAsFactors=F, quote="\"")

colnames(geneinfo)=c("ID", "Name", "Chr", "Start", "End", "Strand", "Biotype")
rownames(geneinfo)=geneinfo$ID

geneinfo=geneinfo[which(geneinfo$Chr%in%c(as.character(1:22), "X", "Y")),]

###############################################################################

lnc=geneinfo$ID[which(geneinfo$Biotype=="lncRNA")]
pc=geneinfo$ID[which(geneinfo$Biotype=="protein_coding")]

###############################################################################

save(list=c("geneinfo", "lnc", "pc"), file=paste(pathRData, "data.gene.info.RData",sep=""))

###############################################################################
