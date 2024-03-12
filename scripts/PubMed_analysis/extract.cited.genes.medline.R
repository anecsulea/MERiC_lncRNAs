###############################################################################

pathPubMed="../../results/PubMed_search/"
pathEnsembl="../../data/ensembl_annotations/"

release=109

###############################################################################

geneinfo=read.table(paste(pathEnsembl, "GeneInfo_Ensembl",release, ".txt", sep=""), h=T, sep="\t", stringsAsFactors=F, quote="\"")

colnames(geneinfo)=c("ID", "Name", "Chr", "Start", "End", "Strand", "Biotype")
rownames(geneinfo)=geneinfo$ID

###############################################################################

ensembl.lnc=geneinfo$ID[which(geneinfo$Biotype=="lncRNA")]
ensembl.pc=geneinfo$ID[which(geneinfo$Biotype=="protein_coding")]

###############################################################################

articles=read.table(paste(pathPubMed, "PubMed_hepatocellular_carcinoma_Title_11_03_2024.tsv", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")

print(paste(nrow(articles), "articles initially"))

###############################################################################
## filter entries

## keep only "Journal Article" and "Retracted Publication"

articles=articles[which(articles$PublicationType%in%c("Journal Article", "Retracted Publication")),]

print(paste(nrow(articles), "articles after filtering types"))

## keep only articles that cite genes

articles=articles[which(articles$CitedGenes!=""),]

print(paste(nrow(articles), "articles after removing those that do not cite genes"))

###############################################################################

## check if article cites lncRNAs

articles$CitedLnc=unlist(lapply(articles$CitedGenes, function(x) paste(intersect(unlist(strsplit(x, split=";")), ensembl.lnc), collapse=";")))
