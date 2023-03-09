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

articles=readLines(paste(pathPubMed, "formatted_results_hepatocellular_carcinoma_Title.txt", sep=""))

###############################################################################

## each article has 8 fields

pmid=unlist(lapply(articles[seq(from=1, to=length(articles), by=8)], function(x) unlist(strsplit(x, split="\\: "))[2]))
hccarticle.genes=unlist(lapply(articles[seq(from=6, to=length(articles), by=8)], function(x) unlist(strsplit(x, split="\\: "))[2]))
cited.genes=unique(unlist(lapply(hccarticle.genes, function(x) {y=unlist(strsplit(x, split=" ")); return(y[seq(from=1, to=length(y), by=2)])})))

cited.pc=intersect(cited.genes, ensembl.pc)
cited.lnc=intersect(cited.genes, ensembl.lnc)

###############################################################################

writeLines(cited.pc, con=paste(pathPubMed, "cited_protein_coding_hepatocellular_carcinoma_Title.txt", sep=""))

writeLines(cited.lnc, con=paste(pathPubMed, "cited_lncRNA_hepatocellular_carcinoma_Title.txt", sep=""))

###############################################################################
