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

###############################################################################

pmid=unlist(lapply(articles[seq(from=1, to=length(articles), by=8)], function(x) unlist(strsplit(x, split="\\: "))[2]))

hccarticle.year=as.numeric(unlist(lapply(articles[seq(from=4, to=length(articles), by=8)], function(x) unlist(strsplit(x, split="\\: "))[2])))


hccarticle.genes=unlist(lapply(articles[seq(from=6, to=length(articles), by=8)], function(x) unlist(strsplit(x, split="\\: "))[2]))
hccarticle.genes=lapply(hccarticle.genes, function(x) {y=unlist(strsplit(x, split=" ")); return(y[seq(from=1, to=length(y), by=2)])})

hccarticle.lnc=unlist(lapply(articles[seq(from=7, to=length(articles), by=8)], function(x) unlist(strsplit(x, split="\\: "))[2]))

cited.genes=unique(unlist(hccarticle.genes))
cited.genes=setdiff(cited.genes, NA)

cited.pc=intersect(cited.genes, ensembl.pc)
cited.lnc=intersect(cited.genes, ensembl.lnc)

###############################################################################

writeLines(cited.pc, con=paste(pathPubMed, "cited_protein_coding_hepatocellular_carcinoma_Title.txt", sep=""))

writeLines(cited.lnc, con=paste(pathPubMed, "cited_lncRNA_hepatocellular_carcinoma_Title.txt", sep=""))

###############################################################################

nb.articles.per.year=as.numeric(table(hccarticle.year))
names(nb.articles.per.year)=as.character(levels(as.factor(hccarticle.year)))

fraction.lncarticles.per.year=tapply(hccarticle.lnc, as.factor(hccarticle.year), function(x) length(which(x=="True"))/length(x))
names(fraction.lncarticles.per.year)=as.character(levels(as.factor(hccarticle.year)))

nb.lncarticles.per.year=tapply(hccarticle.lnc, as.factor(hccarticle.year), function(x) length(which(x=="True")))
names(nb.lncarticles.per.year)=as.character(levels(as.factor(hccarticle.year)))

###############################################################################

all.cited.genes=unlist(hccarticle.genes)
all.cited.genes=all.cited.genes[which(!is.na(all.cited.genes))]

nb.articles.per.gene=as.numeric(table(all.cited.genes))
names(nb.articles.per.gene)=levels(as.factor(all.cited.genes))

###############################################################################

results=data.frame("GeneID"=cited.genes, "GeneName"=geneinfo[cited.genes,"Name"], "GeneType"=geneinfo[cited.genes,"Biotype"], "NbCitations"=nb.articles.per.gene[cited.genes])

results=results[order(results$NbCitations, decreasing=T),]

write.table(results, paste(pathPubMed, "number_of_citations_hepatocellular_carcinoma_Title.txt", sep=""), row.names=F, col.names=T, sep="\t", quote=F)

###############################################################################

## retracted articles

retracted=readLines(paste(pathPubMed, "formatted_results_retracted_hepatocellular_carcinoma_Title.txt", sep=""))

###############################################################################

## each article has 8 fields

pmid.retracted=unlist(lapply(retracted[seq(from=1, to=length(retracted), by=8)], function(x) unlist(strsplit(x, split="\\: "))[2]))

retractedarticle.lnc=unlist(lapply(retracted[seq(from=7, to=length(retracted), by=8)], function(x) unlist(strsplit(x, split="\\: "))[2]))
retractedarticle.year=unlist(lapply(retracted[seq(from=4, to=length(retracted), by=8)], function(x) unlist(strsplit(x, split="\\: "))[2]))

###############################################################################
