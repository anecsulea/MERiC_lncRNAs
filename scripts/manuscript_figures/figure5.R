##########################################################################

if(!("pathFigures"%in%objects())){
    pathRData="../../data_for_publication/RData/"
    pathFigures="../../data_for_publication/figures/"

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

    years=2009:2022

    nb.articles.per.year=table(factor(articles$Year, levels=years))
    nb.retractions.per.year=table(factor(retractions$Year, levels=years))
    nb.pubs.per.year=nb.articles.per.year+nb.retractions.per.year

    nb.lncRNA.articles.per.year=table(factor(articles$Year[which(articles$CitedLnc!="")], levels=years))
    nb.lncRNA.retractions.per.year=table(factor(retractions$Year[which(retractions$CitedLnc!="")], levels=years))
    nb.lncRNA.pubs.per.year=nb.lncRNA.articles.per.year+nb.lncRNA.retractions.per.year

    nb.pc.articles.per.year=table(factor(articles$Year[which(articles$CitedPc!="")], levels=years))
    nb.pc.retractions.per.year=table(factor(retractions$Year[which(retractions$CitedPc!="")], levels=years))
    nb.pc.pubs.per.year=nb.pc.articles.per.year+nb.pc.retractions.per.year

    nb.other.articles.per.year=table(factor(articles$Year[which(articles$CitedLnc=="" & articles$CitedPc=="")], levels=years))
    nb.other.retractions.per.year=table(factor(retractions$Year[which(retractions$CitedLnc=="" & retractions$CitedPc=="")], levels=years))
    nb.other.pubs.per.year=nb.other.articles.per.year+nb.other.retractions.per.year

    ## proportion of lncRNA-citing or pc-citing articles that get retracted

    prop.retracted.lncRNA.articles.per.year=nb.lncRNA.retractions.per.year/(nb.lncRNA.retractions.per.year+nb.lncRNA.articles.per.year)
    prop.retracted.pc.articles.per.year=nb.pc.retractions.per.year/(nb.pc.retractions.per.year+nb.pc.articles.per.year)
    prop.retracted.other.articles.per.year=nb.other.retractions.per.year/(nb.other.retractions.per.year+nb.other.articles.per.year)

    m.retracted=matrix(100*c(prop.retracted.pc.articles.per.year, prop.retracted.lncRNA.articles.per.year, prop.retracted.other.articles.per.year), nrow=3, byrow=T)
    rownames(m.retracted)=c("pc", "lnc", "other")
    colnames(m.retracted)=as.character(years)

    ## are lncRNA-citing articles more often retracted?

    nb.retracted.lnc=sum(nb.lncRNA.retractions.per.year[as.character(years)])
    nb.tot.lnc=sum(nb.lncRNA.retractions.per.year[as.character(years)]+nb.lncRNA.articles.per.year[as.character(years)])

    nb.retracted.pc=sum(nb.pc.retractions.per.year[as.character(years)])
    nb.tot.pc=sum(nb.pc.retractions.per.year[as.character(years)]+nb.pc.articles.per.year[as.character(years)])

    nb.retracted.other=sum(nb.other.retractions.per.year[as.character(years)])
    nb.tot.other=sum(nb.other.retractions.per.year[as.character(years)]+nb.other.articles.per.year[as.character(years)])


    prepare=FALSE
}

###############################################################################

pdf(file=paste(pathFigures, "Figure5.pdf",sep=""), width=6.855, height=3.5)

###############################################################################

## layout

m=matrix(rep(NA, 1*19), nrow=1)

for(i in 1){
    m[i,]=c(rep(1, 4), rep(2, 15))
}

layout(m)

###############################################################################

## proportion retracted articles

par(mar=c(5.5, 3.5, 2, 1))

percent.retracted.bycategory=100*c(nb.retracted.pc, nb.retracted.lnc, nb.retracted.other)/c(nb.tot.pc, nb.tot.lnc, nb.tot.other)

xpos=1:3
width=0.25

xlim=c(0.25, 3.75)
ylim=c(0, max(percent.retracted.bycategory)+1)

plot(1, xlim=xlim, ylim=ylim, type="n", axes=F, xlab="", ylab="")
rect(xpos-width, 0,  xpos+width, percent.retracted.bycategory, col=c("indianred", "steelblue", "gray40"), border=NA)

axis(side=2, mgp=c(3, 0.65, 0))
axis(side=1, at=xpos, label=rep("", length(xpos)), mgp=c(3, 0.5, 0))

mtext("% retracted articles", side=2, line=2.1, cex=0.75)

legend("bottomleft", legend=c("articles citing protein-coding genes"), fill="indianred", bty="n", cex=1.1, horiz=T, inset=c(-0.6, -0.18), xpd=NA)
legend("bottomleft", legend=c("articles citing lncRNAs"), fill="steelblue", bty="n", cex=1.1, horiz=T, inset=c(-0.6, -0.24), xpd=NA)
legend("bottomleft", legend=c("other articles"), fill="gray40", bty="n", cex=1.1, horiz=T, inset=c(-0.6, -0.3), xpd=NA)

mtext("a", side=3, line=0.5, at=-1.45, font=2, cex=0.9)

###############################################################################

par(mar=c(5, 5.1, 2, 1))

b=barplot(m.retracted, beside=T, col=c("indianred", "steelblue", "gray40"), border=c("indianred", "steelblue", "gray40"), names=rep("", ncol(m.retracted)), axes=F)

axis(side=2, mgp=c(3, 0.65, 0))
xpos=apply(b, 2, mean)
axis(side=1, at=xpos, label=rep("", length(xpos)), mgp=c(3, 0.5, 0))

mtext("year of publication", side=1, line=3.5, cex=0.75)
mtext(years, at=xpos, side=1, las=2, line=0.75, cex=0.7)

mtext("% retracted articles", side=2, line=2.25, cex=0.75)

mtext("b", font=2, side=3, at=-5.9, line=0.5, cex=1)

###############################################################################

dev.off()

###############################################################################
