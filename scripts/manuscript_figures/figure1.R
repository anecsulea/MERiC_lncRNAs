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

    ## nb articles per year

    nb.articles.per.year=as.numeric(table(pubs$Year))
    names(nb.articles.per.year)=as.character(levels(as.factor(pubs$Year)))

    ## nb lncRNA articles per year

    nb.lncarticles.per.year=as.numeric(table(pubs$Year[which(pubs$CitedLnc!="")]))
    names(nb.lncarticles.per.year)=as.character(levels(as.factor(pubs$Year[which(pubs$CitedLnc!="")])))

    prepare=FALSE
}

##########################################################################
##########################################################################

pdf(paste(pathFigures, "Figure1.pdf", sep=""), width=6.85, height=3.5)

##########################################################################

## layout

m=matrix(rep(NA, 2*20), nrow=2)

for(i in 1:2){
  m[i,]=c(rep(1, 8), rep(2, 12))
}

layout(m)

##########################################################################

## number of lncRNAs articles per year

years=2009:2022
nbtot=nb.articles.per.year[as.character(years)]
nblnc=nb.lncarticles.per.year[as.character(years)]
proplnc=100*nblnc/nbtot

xpos=1:length(years)

width=0.25

xlim=c(0.5, length(years)+0.5)
ylim=c(0, max(proplnc)*1.1)

par(mar=c(6.0, 3.75, 2.1, 0.75))

plot(1, xlim=xlim, ylim=ylim, type="n", axes=F, xlab="", ylab="")
rect(xpos-width, 0,  xpos+width, proplnc, col="steelblue", border=NA)

## nb lncRNA papers each year

stx=c(1:length(nblnc))
text(nblnc[stx], x=xpos[stx], y=proplnc[stx]+0.35, cex=1)


axis(side=1, at=xpos, labels=rep("", length(xpos)), mgp=c(3, 0.5, 0))
mtext(years, at=xpos, line=0.85,  adj=1, side=1, las=2, cex=0.7)

mtext("year of publication", side=1, line=3.75, cex=0.85)

axis(side=2, cex.axis=1.05, mgp=c(3, 0.75, 0))
mtext("% of all HCC publications", side=2, line=2.25, cex=0.85)

mtext("a", side=3, line=0.5, at=-2.6, font=2, cex=1.1)

#############################################################################

## number of articles per lnc

nbmax=max(nb.citations.lnc)
nb.articles.lnc=table(factor(nb.citations.lnc, levels=1:nbmax))

xpos=1:length(nb.articles.lnc)
width=0.25

xlim=c(0.75, 52.25)
ylim=c(0, 300)

par(mar=c(4.75, 3.5, 2.1, 0.5))

plot(1, xlim=xlim, ylim=ylim, type="n", axes=F, xlab="", ylab="")
rect(xpos-width, 0,  xpos+width, nb.articles.lnc, col="steelblue", border=NA)


xax=c(1,seq(from=5, to=50, by=5))
axis(side=1, at=xax, labels=rep("", length(xax)), mgp=c(3, 0.5, 0), cex.axis=0.95)
mtext(xax, at=xax, line=0.85,  side=1, cex=0.8)
mtext("number of publications citing each lncRNA", side=1, cex=0.85, line=2.25)

axis(side=2, mgp=c(3, 0.75, 0), cex.axis=1.05, at=seq(from=0, to=300, by=50))
mtext("number of lncRNAs", side=2, cex=0.85, line=2.25)

## legend for lncRNA names, for those lncRNAs cited at least 20 times

allnb=nb.citations.lnc
uniquenb=sort(unique(allnb[which(allnb>=20)]))

xtext=uniquenb
ytext=c(85, 145, 115, 85, 55, 115, 55)

names(xtext)=as.character(uniquenb)
names(ytext)=as.character(uniquenb)


## first we draw the arrows
for(i in uniquenb){
  this.lnc=sort(geneinfo[names(nb.citations.lnc)[which(nb.citations.lnc==i)],"Name"])
  nb=length(this.lnc)

  nbpairs=ceiling(nb/2)

  this.ytext=ytext[as.character(i)]

  this.ystartrect=this.ytext-5*nbpairs
  arrows(i, this.ystartrect-2, i, nb+2, length=0.05, lwd=1.25)
}

## then rectangles and text

rectmargin=1

for(i in uniquenb){
  this.lnc=sort(geneinfo[names(nb.citations.lnc)[which(nb.citations.lnc==i)],"Name"])
  nb=length(this.lnc)
  nbpairs=ceiling(nb/2)

  label=paste(unlist(lapply(1:nbpairs, function(x) {y=(2*(x-1)+1):(2*x); y=intersect(y, 1:nb);  return(paste(this.lnc[y], collapse=", "))})), collapse="\n");

  maxchar=max(unlist(lapply(1:nbpairs, function(x) {y=(2*(x-1)+1):(2*x); y=intersect(y, 1:nb);  return(nchar(paste(this.lnc[y], collapse=", ")))})))

  this.xtext=xtext[as.character(i)]
  this.ytext=ytext[as.character(i)]

  this.xstartrect=this.xtext-rectmargin-maxchar/2.35
  this.xendrect=this.xtext+rectmargin+maxchar/2.35

  if(this.xstartrect<1.35){
    this.xstartrect=1.35
    this.xendrect=2*this.xtext-this.xstartrect
  }

  this.ystartrect=this.ytext-9-8*(nbpairs-1)
  this.yendrect=this.ytext+9+8*(nbpairs-1)

  rect(this.xstartrect, this.ystartrect, this.xendrect, this.yendrect, xpd=NA, col="white")

  text(label, x=this.xtext, y=this.ytext, cex=1)

}

mtext("b", side=3, line=0.5, at=-7.15, font=2, cex=1.1)

##########################################################################

dev.off()

##########################################################################
##########################################################################


