###############################################################################

pathPubMed="../../results/PubMed_search/"
pathEnsembl="../../data/ensembl_annotations/"

release=109

###############################################################################

if(redo){
geneinfo=read.table(paste(pathEnsembl, "GeneInfo_Ensembl",release, ".txt", sep=""), h=T, sep="\t", stringsAsFactors=F, quote="\"")

colnames(geneinfo)=c("ID", "Name", "Chr", "Start", "End", "Strand", "Biotype")
rownames(geneinfo)=geneinfo$ID

###############################################################################

ensembl.lnc=geneinfo$ID[which(geneinfo$Biotype=="lncRNA")]
ensembl.pc=geneinfo$ID[which(geneinfo$Biotype=="protein_coding")]

###############################################################################

pubs=read.table(paste(pathPubMed, "PubMed_hepatocellular_carcinoma_Title_11_03_2024.tsv", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")

print(paste(nrow(pubs), "publications initially"))

###############################################################################
## filter entries

## keep only "Journal Article" and "Retracted Publication"

pubs=pubs[which(pubs$PublicationType%in%c("Journal Article", "Retracted Publication")),]

print(paste(nrow(pubs), "pubs after filtering types"))

## keep only pubs that cite genes

pubs=pubs[which(pubs$CitedGenes!=""),]

print(paste(nrow(pubs), "pubs after removing those that do not cite genes"))

###############################################################################

## keep only pubs between 2000 and 2022
## 2023 : very few retractions, they may not have had time

pubs$Year=as.numeric(pubs$Year)
pubs=pubs[which(pubs$Year>=2000 & pubs$Year<=2022),]

###############################################################################

## check if article cites lncRNAs

pubs$CitedLnc=unlist(lapply(pubs$CitedGenes, function(x) paste(intersect(unlist(strsplit(x, split=";")), ensembl.lnc), collapse=";")))
pubs$CitedPc=unlist(lapply(pubs$CitedGenes, function(x) paste(intersect(unlist(strsplit(x, split=";")), ensembl.pc), collapse=";")))

###############################################################################

## divide into articles and retracted articles

articles=pubs[which(pubs$PublicationType=="Journal Article"),]
retractions=pubs[which(pubs$PublicationType=="Retracted Publication"),]

###############################################################################

## statistics per year

years=2000:2022

nb.articles.per.year=table(factor(articles$Year, levels=2000:2022))
nb.retractions.per.year=table(factor(retractions$Year, levels=2000:2022))
nb.pubs.per.year=nb.articles.per.year+nb.retractions.per.year

nb.other.articles.per.year=table(factor(articles$Year[which(articles$CitedLnc=="")], levels=2000:2022))
nb.other.retractions.per.year=table(factor(retractions$Year[which(retractions$CitedLnc=="")], levels=2000:2022))
nb.other.pubs.per.year=nb.other.articles.per.year+nb.other.retractions.per.year

nb.lncRNA.articles.per.year=table(factor(articles$Year[which(articles$CitedLnc!="")], levels=2000:2022))
nb.lncRNA.retractions.per.year=table(factor(retractions$Year[which(retractions$CitedLnc!="")], levels=2000:2022))
nb.lncRNA.pubs.per.year=nb.lncRNA.articles.per.year+nb.lncRNA.retractions.per.year

## proportion of publications that cite lncRNA

prop.pubs.lncRNAs.per.year=nb.lncRNA.pubs.per.year/nb.pubs.per.year

## proportion of not-retracted articles that cite lncRNAs
prop.articles.lncRNAs.per.year=nb.lncRNA.articles.per.year/nb.articles.per.year

## proportion of retracted articles that cite lncRNAs
prop.retractions.lncRNAs.per.year=nb.lncRNA.retractions.per.year/nb.retractions.per.year

## proportion of lncRNA-citing articles that get retracted

prop.retracted.lncRNA.articles.per.year=nb.lncRNA.retractions.per.year/(nb.lncRNA.retractions.per.year+nb.lncRNA.articles.per.year)

prop.retracted.other.articles.per.year=nb.other.retractions.per.year/(nb.other.retractions.per.year+nb.other.articles.per.year)

prop.retracted.pubs.per.year=nb.retractions.per.year/nb.pubs.per.year

}
###############################################################################

## are lncRNA-citing articles more often retracted?

nb.retracted.lnc=sum(nb.lncRNA.retractions.per.year[2000:2022])
nb.tot.lnc=sum(nb.lncRNA.retractions.per.year+nb.lncRNA.articles.per.year)

nb.retracted.other=sum(nb.other.retractions.per.year)
nb.tot.other=sum(nb.other.retractions.per.year+nb.other.articles.per.year)

###############################################################################

pdf(file="Figure1.pdf", width=8.5, height=6.5)

###############################################################################

par(mfrow=c(2,2))

par(mar=c(4,3.5,2.5,3.5))

###############################################################################

## nb pubs per year

ylim=c(0,3500)

b=barplot(nb.pubs.per.year, col="lightblue",axes=F, names=rep("", length(years)), ylim=ylim, space=0.5, xlim=c(0,37))

## add total number
nb.tot.pubs=sum(nb.pubs.per.year)

ratio=max(nb.pubs.per.year)/nb.tot.pubs

xpos=36.5
width=0.5

rect(xpos-width,0, xpos+width, nb.tot.pubs*ratio, col="gray40")

yax=pretty(c(0, nb.tot.pubs))
yax=yax[which(yax<(ylim[2]/ratio))]

axis(side=4, at=yax*ratio, labels=yax, mgp=c(3,0.65,0), cex.axis=0.9)

## normal axes

axis(side=2, mgp=c(3,0.65,0), cex.axis=0.9)
axis(side=1, at=c(b, 36.5), labels=rep("", length(b)+1))

mtext(years, side=1, at=b, line=0.75, las=2, cex=0.7)
mtext("total", side=1, at=36.5, line=0.8, las=2, cex=0.72)

mtext("year of publication", side=1, line=2.5, cex=0.8)
mtext("number of publications", side=2, line=2.25, cex=0.8)

mtext("a", font=2, side=3, at=-7.5, line=1.5, cex=1)

###############################################################################

## % of papers that cite lncRNAs

b=barplot(100*prop.pubs.lncRNAs.per.year,col="lightblue",axes=F, names=rep("", length(years)), space=0.5, xlim=c(0,37))

## add total proportion
prop.tot=100*sum(nb.lncRNA.pubs.per.year)/sum(nb.pubs.per.year)

xpos=36.5
width=0.5

rect(xpos-width,0, xpos+width, prop.tot, col="gray40")

axis(side=2, mgp=c(3,0.65,0), cex.axis=0.9)
axis(side=1, at=c(b,36.5), labels=rep("", length(b)+1))

mtext(years, side=1, at=b, line=0.75, las=2, cex=0.7)
mtext("total", side=1, at=36.5, line=0.8, las=2, cex=0.72)

mtext("year of publication", side=1, line=2.5, cex=0.8)
mtext("% publications that cite lncRNAs", side=2, line=2.25, cex=0.8)

mtext("b", font=2, side=3, at=-7.5, line=1.5, cex=1)

###############################################################################

## publications that are retracted every year

b=barplot(100*prop.retracted.pubs.per.year,col="lightblue",axes=F, names=rep("", length(years)), ylim=c(0,6.5), space=0.5, xlim=c(0,37))

## add total proportion
prop.tot=100*sum(nb.retractions.per.year)/sum(nb.pubs.per.year)

xpos=36.5
width=0.5

rect(xpos-width,0, xpos+width, prop.tot, col="gray40")

axis(side=2, mgp=c(3,0.65,0), cex.axis=0.9)
axis(side=1, at=c(b,36.5), labels=rep("", length(b)+1))
mtext(years, side=1, at=b, line=0.75, las=2, cex=0.7)

mtext("year of publication", side=1, line=2.5, cex=0.8)
mtext("total", side=1, at=36.5, line=0.8, las=2, cex=0.72)

mtext("% retracted publications", side=2, line=2.25, cex=0.8)

mtext("c", font=2, side=3, at=-7.5, line=1.5, cex=1)

###############################################################################

## lncRNA-citing publications that are retracted every year

b=barplot(100*prop.retracted.lncRNA.articles.per.year,col="lightblue",axes=F, names=rep("", length(years)), ylim=c(0,6.5), space=0.5, xlim=c(0,37))


## add total proportion
prop.tot=100*sum(nb.lncRNA.retractions.per.year)/sum(nb.lncRNA.pubs.per.year)

xpos=36.5
width=0.5

rect(xpos-width,0, xpos+width, prop.tot, col="gray40")


axis(side=2, mgp=c(3,0.65,0), cex.axis=0.9)
axis(side=1, at=c(b, 36.5), labels=rep("", length(b)+1))

mtext(years, side=1, at=b, line=0.75, las=2, cex=0.7)
mtext("total", side=1, at=36.5, line=0.8, las=2, cex=0.72)

mtext("year of publication", side=1, line=2.5, cex=0.8)
mtext("% retracted lncRNA publications", side=2, line=2.25, cex=0.8)

mtext("d", font=2, side=3, at=-7.5, line=1.5, cex=1)

###############################################################################

dev.off()

###############################################################################
###############################################################################

## for each lncRNA, how many times was it cited

all.cited.lnc=unlist(lapply(pubs$CitedLnc, function(x) setdiff(unlist(strsplit(x, split=";")), "")))
nb.citations.per.lnc=table(all.cited.lnc)

retracted.cited.lnc=unlist(lapply(pubs$CitedLnc[which(pubs$PublicationType=="Retracted Publication")], function(x) setdiff(unlist(strsplit(x, split=";")), "")))
nb.retracted.citations.per.lnc=table(retracted.cited.lnc)

###############################################################################
