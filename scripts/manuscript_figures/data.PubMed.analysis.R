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

###############################################################################

## check if article cites lncRNAs or protein-coding genes

pubs$CitedLnc=unlist(lapply(pubs$CitedGenes, function(x) paste(intersect(unlist(strsplit(x, split=";")), lnc), collapse=";")))
pubs$CitedPc=unlist(lapply(pubs$CitedGenes, function(x) paste(intersect(unlist(strsplit(x, split=";")), pc), collapse=";")))

##########################################################################


##########################################################################



