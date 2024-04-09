##########################################################################

if(!("pathFigures"%in%objects())){
    pathRData="../../data_for_publication/RData/"
    pathFigures="../../data_for_publication/main_figures/"

    load=TRUE
    prepare=TRUE
}

##########################################################################

if(load){
    load(paste(pathRData, "data.PubMed.analysis.RData", sep=""))

    load(paste(pathRData, "data.gene.info.RData", sep=""))

    load=FALSE
}

##########################################################################

if(prepare){
    ## analyze publications

    articles=pubs[which(pubs$PublicationType=="Journal Article"),]
    retractions=pubs[which(pubs$PublicationType=="Retracted Publication"),]

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

    ## are lncRNA-citing articles more often retracted?

    nb.retracted.lnc=sum(nb.lncRNA.retractions.per.year[2000:2022])
    nb.tot.lnc=sum(nb.lncRNA.retractions.per.year+nb.lncRNA.articles.per.year)

    nb.retracted.other=sum(nb.other.retractions.per.year)
    nb.tot.other=sum(nb.other.retractions.per.year+nb.other.articles.per.year)

## stats for retractions

    years=2009:2022

    prop.retracted=100*prop.retracted.pubs.per.year
    prop.retracted=prop.retracted[as.character(years)]

    prop.tot.retracted=100*sum(nb.retractions.per.year[as.character(years)])/sum(nb.pubs.per.year[as.character(years)])


    prop.retracted.lnc=100*prop.retracted.lncRNA.articles.per.year
    prop.retracted.lnc=prop.retracted.lnc[as.character(years)]

    ## add total proportion
    prop.tot.retracted.lnc=100*sum(nb.lncRNA.retractions.per.year[as.character(years)])/sum(nb.lncRNA.pubs.per.year[as.character(years)])


    prepare=FALSE
}

###############################################################################

pdf(file=paste(pathFigures, "Figure5.pdf",sep=""), width=6.855, height=3.5)

###############################################################################

## layout

m=matrix(rep(NA, 1*10), nrow=1)

for(i in 1){
    m[i,]=c(rep(1, 5), rep(2, 5))
}

layout(m)


###############################################################################

par(mar=c(5,3.5,2.5,1.5))

years=2009:2022

xpos=1:length(years)

width=0.25

xlim=c(0.5, length(years)+2)
ylim=c(0, max(c(prop.retracted, prop.retracted.lnc))*1.1)


###############################################################################

## publications that are retracted every year

plot(1, xlim=xlim, ylim=ylim, type="n", axes=F, xlab="", ylab="")
rect(xpos-width, 0,  xpos+width, prop.retracted, col="steelblue", border=NA)

## add total proportion

rect(max(xpos)+1.5-width,0, max(xpos)+1.5+width, prop.tot.retracted, col="gray40", border=NA)

axis(side=2, mgp=c(3,0.65,0), cex.axis=0.95)
axis(side=1, at=c(xpos, max(xpos)+1.5), labels=rep("", length(xpos)+1))
mtext(years, side=1, at=xpos, line=0.75, las=2, cex=0.7)

mtext("year of publication", side=1, line=3.5, cex=0.8)
mtext("total", side=1, at=max(xpos)+1.5, line=0.8, las=2, cex=0.72)

mtext("% retracted publications", side=2, line=2.25, cex=0.8)

mtext("a", font=2, side=3, at=-2.25, line=0.5, cex=1)

###############################################################################

## lncRNA-citing publications that are retracted every year

plot(1, xlim=xlim, ylim=ylim, type="n", axes=F, xlab="", ylab="")
rect(xpos-width, 0,  xpos+width, prop.retracted.lnc, col="steelblue", border=NA)


rect(max(xpos)+1.5-width,0, max(xpos)+1.5+width, prop.tot.retracted.lnc, col="gray40", border=NA)

axis(side=2, mgp=c(3,0.65,0), cex.axis=0.95)
axis(side=1, at=c(xpos, max(xpos)+1.5), labels=rep("", length(xpos)+1))
mtext(years, side=1, at=xpos, line=0.75, las=2, cex=0.7)

mtext("year of publication", side=1, line=3.5, cex=0.8)
mtext("total", side=1, at=max(xpos)+1.5, line=0.8, las=2, cex=0.72)

mtext("% retracted lncRNA publications", side=2, line=2.25, cex=0.8)

mtext("b", font=2, side=3, at=-2.25, line=0.5, cex=1)

###############################################################################

dev.off()

###############################################################################
