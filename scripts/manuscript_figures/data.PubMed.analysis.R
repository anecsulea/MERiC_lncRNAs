##########################################################################

pathPubMed="../../results/PubMed_search/"
pathRData="../../data_for_publication/RData/"

##########################################################################

## load gene info
load(paste(pathRData,"data.gene.info.RData",sep=""))

##########################################################################

## all publications
pubs=read.table(paste(pathPubMed, "PubMed_hepatocellular_carcinoma_Title_11_03_2024.tsv", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")

print(paste(nrow(pubs), "publications initially"))

## keep only journal articles and retracted articles

pubs=pubs[which(pubs$PublicationType%in%c("Journal Article", "Retracted Publication")),]

print(paste(nrow(pubs), "pubs after filtering types"))

## keep only pubs that cite genes

pubs=pubs[which(pubs$CitedGenes!=""),]

print(paste(nrow(pubs), "pubs after removing those that do not cite genes"))

pubs$Year=as.numeric(pubs$Year)
pubs=pubs[which(pubs$Year>=2000 & pubs$Year<=2022),]

print(paste(nrow(pubs), "pubs after filtering years"))

###############################################################################

## check if article cites lncRNAs or protein-coding genes

pubs$CitedLnc=unlist(lapply(pubs$CitedGenes, function(x) paste(intersect(unlist(strsplit(x, split=";")), lnc), collapse=";")))
pubs$CitedPc=unlist(lapply(pubs$CitedGenes, function(x) paste(intersect(unlist(strsplit(x, split=";")), pc), collapse=";")))

##########################################################################

## extract all cited pc and lnc

all.cited.pc=unlist(lapply(pubs$CitedPc, function(x) unlist(strsplit(x, split=","))))
nb.citations.pc=as.numeric(table(all.cited.pc))
names(nb.citations.pc)=levels(as.factor(all.cited.pc))

all.cited.lnc=unlist(lapply(pubs$CitedLnc, function(x) unlist(strsplit(x, split=","))))
nb.citations.lnc=as.numeric(table(all.cited.lnc))
names(nb.citations.lnc)=levels(as.factor(all.cited.lnc))

##########################################################################

save(list=c("pubs", "all.cited.pc", "all.cited.lnc", "nb.citations.pc", "nb.citations.lnc"), file=paste(pathRData, "data.PubMed.analysis.RData",sep=""))

##########################################################################


