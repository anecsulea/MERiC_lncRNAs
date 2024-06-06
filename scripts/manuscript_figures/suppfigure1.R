##########################################################################

if(!("pathFigures"%in%objects())){
    pathRData="../../data_for_publication/RData/"
    pathFigures="../../data_for_publication/figures/"

    load=TRUE
    prepare=TRUE
}

##########################################################################

if(load){
    load(paste(pathRData, "data.PubMed.analysis.2023.RData", sep=""))

    load(paste(pathRData, "data.gene.info.RData", sep=""))

    load=FALSE
}

##########################################################################

if(prepare){

    ## nb articles per year

    nb.articles.per.year=as.numeric(table(pubs$Year))
    names(nb.articles.per.year)=as.character(levels(as.factor(pubs$Year)))

    ## nb protein-coding articles per year

    nb.pcarticles.per.year=as.numeric(table(pubs$Year[which(pubs$CitedPc!="")]))
    names(nb.pcarticles.per.year)=as.character(levels(as.factor(pubs$Year[which(pubs$CitedPc!="")])))

    prepare=FALSE
}

##########################################################################

pdf(paste(pathFigures, "SupplementaryFigure1.pdf", sep=""), width=6.85, height=4)

##########################################################################

## layout

m=matrix(rep(NA, 2*20), nrow=2)

for(i in 1:2){
  m[i,]=c(rep(1, 7), rep(2, 13))
}

layout(m)

##########################################################################

## number of pc articles per year

years=2009:2023
nbtot=nb.articles.per.year[as.character(years)]
nbpc=nb.pcarticles.per.year[as.character(years)]
proppc=100*nbpc/nbtot

xpos=1:length(years)

width=0.25

xlim=c(0.5, length(years)+0.5)
ylim=c(0, max(proppc)*1.1)

par(mar=c(6.0, 3.75, 2.1, 0.75))

plot(1, xlim=xlim, ylim=ylim, type="n", axes=F, xlab="", ylab="")
rect(xpos-width, 0,  xpos+width, proppc, col="steelblue", border=NA)

axis(side=1, at=xpos, labels=rep("", length(xpos)), mgp=c(3, 0.5, 0))
mtext(years, at=xpos, line=0.85,  adj=1, side=1, las=2, cex=0.7)

mtext("year of publication", side=1, line=3.75, cex=0.85)

axis(side=2, cex.axis=1.05, mgp=c(3, 0.75, 0))
mtext("% of all HCC publications that cite p-c genes", side=2, line=2.25, cex=0.85)

mtext("a", side=3, line=0.5, at=-3.1, font=2, cex=1.1)

#############################################################################

## number of articles per pc

nbmax=max(nb.citations.pc)
nb.articles.pc=table(cut(nb.citations.pc, breaks=c(0:49,nbmax)))

xpos=1:length(nb.articles.pc)
width=0.25

xlim=c(0.75, 52.25)
ylim=c(0, 3000)

par(mar=c(4.75, 3.5, 2.1, 0.5))

plot(1, xlim=xlim, ylim=ylim, type="n", axes=F, xlab="", ylab="")
rect(xpos-width, 0,  xpos+width, nb.articles.pc, col="steelblue", border=NA)


xax=c(1,seq(from=5, to=50, by=5))
axis(side=1, at=xax, labels=rep("", length(xax)), mgp=c(3, 0.5, 0), cex.axis=0.95)
mtext(xax[-length(xax)], at=xax[-length(xax)], line=0.85,  side=1, cex=0.8)
mtext("50+", at=xax[length(xax)], line=0.85,  side=1, cex=0.8)

mtext("number of publications citing each protein-coding gene", side=1, cex=0.85, line=2.25)

axis(side=2, mgp=c(3, 0.75, 0), cex.axis=1.05, at=seq(from=0, to=3000, by=500))
mtext("number of protein-coding genes", side=2, cex=0.85, line=2.25)

mtext("b", side=3, line=0.5, at=-6.7, font=2, cex=1.1)

##########################################################################

dev.off()

##########################################################################
##########################################################################


