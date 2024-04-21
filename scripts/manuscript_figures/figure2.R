##########################################################################

if(!("pathFigures"%in%objects())){
    pathRData="../../data_for_publication/RData/"
    pathFigures="../../data_for_publication/figures/"

    library(vioplot)

    load=TRUE
    prepare=TRUE
}

##########################################################################

if(load){
    ## cited lncRNAs
    load(paste(pathRData, "data.PubMed.analysis.RData", sep=""))

    ## gene info
    load(paste(pathRData, "data.gene.info.RData", sep=""))

    ## gene overlaps
    load(paste(pathRData, "data.gene.overlaps.RData", sep=""))

    ## TPM
    load(paste(pathRData, "data.expression.levels.MERiC.RData",sep=""))
    tpm.meric=tpm

    ## sample info
    load(paste(pathRData, "data.sample.info.MERiC.RData",sep=""))

    ## sequence conservation
    load(paste(pathRData, "data.phastcons.30way.RData",sep=""))

    load=FALSE
}

##########################################################################

if(prepare){
    lnc.cited.once=names(nb.citations.lnc)[which(nb.citations.lnc==1)]
    lnc.cited.more=names(nb.citations.lnc)[which(nb.citations.lnc>1)]
    lnc.cited.all=names(nb.citations.lnc)
    other.lnc=setdiff(lnc, lnc.cited.all)

    pc.cited.once=names(nb.citations.pc)[which(nb.citations.pc==1)]
    pc.cited.more=names(nb.citations.pc)[which(nb.citations.pc>1)]
    pc.cited.all=names(nb.citations.pc)
    other.pc=setdiff(pc, pc.cited.all)

    ## sense and antisense overlaps
    sense.overlaps=gene.overlaps[which(gene.overlaps$Type=="sense"),]
    antisense.overlaps=gene.overlaps[which(gene.overlaps$Type=="antisense"),]

    ## actual bidirectional promoters
    biprom1kb=biprom1kb[which(!is.na(biprom1kb$GenesCloseTSS)),]
    biprom5kb=biprom5kb[which(!is.na(biprom5kb$GenesCloseTSS)),]

    ########################################################

    prop.antisense.overlaps=list()
    antisense.conf=list()
    prop.biprom=list()
    biprom.conf=list()

    for(genetype in c("lnc.cited.once", "lnc.cited.more", "other.lnc", "pc.cited.once", "pc.cited.more", "other.pc")){
        genes=get(genetype)

        prop.antisense.overlaps[[genetype]]=length(which(genes%in%antisense.overlaps$GeneID))/length(genes)
        antisense.conf[[genetype]]=prop.test(length(which(genes%in%antisense.overlaps$GeneID)), length(genes))$conf
        prop.biprom[[genetype]]=length(which(genes%in%biprom1kb$GeneID))/length(genes)
        biprom.conf[[genetype]]=prop.test(length(which(genes%in%biprom1kb$GeneID)), length(genes))$conf
    }

    ## average expression level across various samples
    meantpm.nontumor=apply(tpm.meric[ ,nontumor.samples$biopsyID],1, mean)
    meantpm.tumor=apply(tpm.meric[ ,tumor.samples$tumor_biopsyID],1, mean)

    prepare=FALSE
}

##########################################################################

pdf(paste(pathFigures, "Figure2.pdf", sep=""), width=6.85, height=5.5)

##########################################################################

## layout

m=matrix(rep(NA, 21*38), nrow=21)

for(i in 1:10){
    m[i,]=c(rep(1, 19), rep(2,19))
}

for(i in 11:20){
    m[i,]=c(rep(3, 12), rep(4, 12), rep(5,14))
}

for(i in 21){
    m[i,]=c(rep(6,38))
}

layout(m)

##########################################################################

genetypes=c("pc.cited.more", "pc.cited.once", "other.pc", "lnc.cited.more", "lnc.cited.once",  "other.lnc")

xpos.genetypes=c(1, 3.5, 6, 2,  4.5, 7)
names(xpos.genetypes)=genetypes

col.genetypes=rep(c("indianred", "steelblue"), each=3)
names(col.genetypes)=genetypes

##########################################################################

## expression levels in meric, tumor samples

par(mar=c(2.5, 3.75, 2.5, 1.75))

xlim=c(0.5,8)
ylim=c(0,20)

plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=ylim)

## expression in tumor samples

for(type in genetypes){
    this.genes=get(type)
    this.xpos=xpos.genetypes[type]
    this.col=col.genetypes[type]

    vioplot(log2(meantpm.tumor[this.genes]+1), h=1, add=T, axes=F, at=this.xpos, border="black", pchMed=21, colMed="black", colMed2="white", cex=0.95, col=this.col)
}

abline(v=mean(xpos.genetypes[c("lnc.cited.more", "pc.cited.once")]), lty=3)
abline(v=mean(xpos.genetypes[c("lnc.cited.once", "other.pc")]), lty=3)

mtext("tumor", side=3, line=0.5, at=mean(xpos.genetypes), cex=0.75)

axis(side=2, mgp=c(3,0.65,0))
mtext("mean expression level (log2 TPM)", side=2, line=2.5, cex=0.75)

mtext("cited", side=1, at=mean(xpos.genetypes[c("pc.cited.more", "lnc.cited.more")]),line=0.5, cex=0.75)
mtext(">1 articles", side=1, at=mean(xpos.genetypes[c("pc.cited.more", "lnc.cited.more")]),line=1.5, cex=0.75)

mtext("cited", side=1, at=mean(xpos.genetypes[c("pc.cited.once", "lnc.cited.once")]),line=0.5, cex=0.75)
mtext("1 article", side=1, at=mean(xpos.genetypes[c("pc.cited.once", "lnc.cited.once")]),line=1.5, cex=0.75)

mtext("not cited", side=1, at=mean(xpos.genetypes[c("other.pc", "other.lnc")]),line=1, cex=0.75)

axis(side=1, at=xpos.genetypes, labels=rep("", length(xpos.genetypes)), mgp=c(3,0.5,0))

legend("topright", legend=c("protein-coding", "lncRNAs"), fill=c("indianred", "steelblue"), xpd=NA, inset=c(-0.05, -0.05), cex=1.1, bty="n")

mtext("a", font=2, line=0.95, at=-1)

##########################################################################

## expression levels in meric, adjacent tissue samples

par(mar=c(2.5, 4.75, 2.5, 0.75))

xlim=c(0.5,8)
ylim=c(0,20)

plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=ylim)

## expression in tumor samples

for(type in genetypes){
    this.genes=get(type)
    this.xpos=xpos.genetypes[type]
    this.col=col.genetypes[type]

    vioplot(log2(meantpm.nontumor[this.genes]+1), h=1, add=T, axes=F, at=this.xpos, border="black", pchMed=21, colMed="black", colMed2="white", cex=0.95, col=this.col)
}

abline(v=mean(xpos.genetypes[c("lnc.cited.more", "pc.cited.once")]), lty=3)
abline(v=mean(xpos.genetypes[c("lnc.cited.once", "other.pc")]), lty=3)

mtext("adjacent tissue", side=3, line=0.5, at=mean(xpos.genetypes), cex=0.75)

axis(side=2, mgp=c(3,0.65,0))
mtext("mean expression level (log2 TPM)", side=2, line=2.5, cex=0.75)

mtext("cited", side=1, at=mean(xpos.genetypes[c("pc.cited.more", "lnc.cited.more")]),line=0.5, cex=0.75)
mtext(">1 articles", side=1, at=mean(xpos.genetypes[c("pc.cited.more", "lnc.cited.more")]),line=1.5, cex=0.75)

mtext("cited", side=1, at=mean(xpos.genetypes[c("pc.cited.once", "lnc.cited.once")]),line=0.5, cex=0.75)
mtext("1 article", side=1, at=mean(xpos.genetypes[c("pc.cited.once", "lnc.cited.once")]),line=1.5, cex=0.75)

mtext("not cited", side=1, at=mean(xpos.genetypes[c("other.pc", "other.lnc")]),line=1, cex=0.75)

axis(side=1, at=xpos.genetypes, labels=rep("", length(xpos.genetypes)), mgp=c(3,0.5,0))

mtext("b", font=2, line=0.95, at=-1)

##########################################################################

## antisense overlap with pc genes

par(mar=c(2.5, 3.75, 3, 1.75))

plot(1, type="n", xlab="", ylab="", axes=F, ylim=c(0,100), xlim=c(0.5,8))

width=diff(xpos.genetypes)[1]/8

for(type in genetypes){
    this.prop=100*prop.antisense.overlaps[[type]]
    this.xpos=xpos.genetypes[type]
    this.col=col.genetypes[type]

    this.conf=antisense.conf[[type]]

    rect(this.xpos-width, 0, this.xpos+width, this.prop, col=this.col)
    segments(this.xpos, 100*this.conf[1], this.xpos, 100*this.conf[2], lwd=1.5)
}

mtext("cited", side=1, at=mean(xpos.genetypes[c("pc.cited.more", "lnc.cited.more")]),line=0.75, cex=0.75)
mtext(">1", side=1, at=mean(xpos.genetypes[c("pc.cited.more", "lnc.cited.more")]),line=1.75, cex=0.75)

mtext("cited", side=1, at=mean(xpos.genetypes[c("pc.cited.once", "lnc.cited.once")]),line=0.75, cex=0.75)
mtext("1", side=1, at=mean(xpos.genetypes[c("pc.cited.once", "lnc.cited.once")]),line=1.75, cex=0.75)

mtext("not cited", side=1, at=mean(xpos.genetypes[c("other.pc", "other.lnc")]),line=1.25, cex=0.75)

axis(side=2, mgp=c(3,0.65,0))

axis(side=1, at=xpos.genetypes, labels=rep("",length(xpos.genetypes)))

abline(v=mean(xpos.genetypes[c("lnc.cited.more", "pc.cited.once")]), lty=3)
abline(v=mean(xpos.genetypes[c("lnc.cited.once", "other.pc")]), lty=3)

mtext("% with antisense overlap", side=2, line=2.5, cex=0.8)

mtext("c", font=2, line=0.95, at=-2)

##########################################################################

## bidirectional promoters

par(mar=c(2.5, 3.25, 3, 2.25))

plot(1, type="n", xlab="", ylab="", axes=F, ylim=c(0,100), xlim=c(0.5,8))

width=diff(xpos.genetypes)[1]/8

for(type in genetypes){
    this.prop=100*prop.biprom[[type]]
    this.xpos=xpos.genetypes[type]
    this.col=col.genetypes[type]

    this.conf=biprom.conf[[type]]

    rect(this.xpos-width, 0, this.xpos+width, this.prop, col=this.col)
    segments(this.xpos, 100*this.conf[1], this.xpos, 100*this.conf[2], lwd=1.5)
}

mtext("cited", side=1, at=mean(xpos.genetypes[c("pc.cited.more", "lnc.cited.more")]),line=0.75, cex=0.75)
mtext(">1", side=1, at=mean(xpos.genetypes[c("pc.cited.more", "lnc.cited.more")]),line=1.75, cex=0.75)

mtext("cited", side=1, at=mean(xpos.genetypes[c("pc.cited.once", "lnc.cited.once")]),line=0.75, cex=0.75)
mtext("1", side=1, at=mean(xpos.genetypes[c("pc.cited.once", "lnc.cited.once")]),line=1.75, cex=0.75)

mtext("not cited", side=1, at=mean(xpos.genetypes[c("other.pc", "other.lnc")]),line=1.25, cex=0.75)

axis(side=2, mgp=c(3,0.65,0))

axis(side=1, at=xpos.genetypes, labels=rep("",length(xpos.genetypes)))

abline(v=mean(xpos.genetypes[c("lnc.cited.more", "pc.cited.once")]), lty=3)
abline(v=mean(xpos.genetypes[c("lnc.cited.once", "other.pc")]), lty=3)

mtext("% with bidirectional promoters", side=2, line=2.5, cex=0.8)

mtext("d", font=2, line=0.95, at=-2)


##########################################################################

## sequence conservation scores

par(mar=c(2.5, 3.25, 3, 0.75))

xlim=c(0.5,8)
ylim=c(0,1)

plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=ylim)

for(type in genetypes){
    this.genes=get(type)
    this.phast=phastcons[this.genes, "Score"]
    this.col=col.genetypes[type]
    this.xpos=xpos.genetypes[type]

    vioplot(this.phast, h=0.02, add=T, axes=F, at=this.xpos, border="black", pchMed=21, colMed="black", colMed2="white", cex=0.95, col=this.col)
}

axis(side=2, mgp=c(3,0.65,0))


axis(side=1, at=xpos.genetypes, labels=rep("",length(xpos.genetypes)))
mtext("sequence conservation", side=2, line=2.5, cex=0.8)



mtext("cited", side=1, at=mean(xpos.genetypes[c("pc.cited.more", "lnc.cited.more")]),line=0.75, cex=0.75)
mtext(">1", side=1, at=mean(xpos.genetypes[c("pc.cited.more", "lnc.cited.more")]),line=1.75, cex=0.75)

mtext("cited", side=1, at=mean(xpos.genetypes[c("pc.cited.once", "lnc.cited.once")]),line=0.75, cex=0.75)
mtext("1", side=1, at=mean(xpos.genetypes[c("pc.cited.once", "lnc.cited.once")]),line=1.75, cex=0.75)

mtext("not cited", side=1, at=mean(xpos.genetypes[c("other.pc", "other.lnc")]),line=1.25, cex=0.75)

mtext("e", font=2, line=0.95, at=-1.35)

##########################################################################

dev.off()

##########################################################################
