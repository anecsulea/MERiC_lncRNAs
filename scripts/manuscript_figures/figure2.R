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
    load(paste(pathRData, "data.PubMed.analysis.2023.RData", sep=""))

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
    lnc.other=setdiff(lnc, lnc.cited.all)

    pc.cited.once=names(nb.citations.pc)[which(nb.citations.pc==1)]
    pc.cited.more=names(nb.citations.pc)[which(nb.citations.pc>1)]
    pc.cited.all=names(nb.citations.pc)
    pc.other=setdiff(pc, pc.cited.all)

    ## sense and antisense overlaps
    sense.overlaps=gene.overlaps[which(gene.overlaps$Type=="sense"),]
    antisense.overlaps=gene.overlaps[which(gene.overlaps$Type=="antisense"),]

    ## actual bidirectional promoters
    biprom1kb=biprom1kb[which(!is.na(biprom1kb$GenesCloseTSS)),]
    biprom5kb=biprom5kb[which(!is.na(biprom5kb$GenesCloseTSS)),]

    ########################################################

    prop.antisense.overlaps=list()
    antisense.overlaps.conf=list()
    prop.biprom=list()
    biprom.conf=list()

    nb.antisense.overlaps=list()
    nb.biprom=list()

    for(genetype in c("lnc.cited.once", "lnc.cited.more", "lnc.other", "pc.cited.once", "pc.cited.more", "pc.other")){
        genes=get(genetype)

        nb.antisense.overlaps[[genetype]]=length(which(genes%in%antisense.overlaps$GeneID))
        nb.biprom[[genetype]]=length(which(genes%in%biprom1kb$GeneID))

        prop.antisense.overlaps[[genetype]]=length(which(genes%in%antisense.overlaps$GeneID))/length(genes)
        antisense.overlaps.conf[[genetype]]=prop.test(length(which(genes%in%antisense.overlaps$GeneID)), length(genes))$conf
        prop.biprom[[genetype]]=length(which(genes%in%biprom1kb$GeneID))/length(genes)
        biprom.conf[[genetype]]=prop.test(length(which(genes%in%biprom1kb$GeneID)), length(genes))$conf
    }

    ## average expression level across various samples
    meantpm.nontumor=apply(tpm.meric[ ,nontumor.samples$biopsyID],1, mean)
    meantpm.tumor=apply(tpm.meric[ ,tumor.samples$tumor_biopsyID],1, mean)

    ## pairwise comparisons

    comparisons=list(c("cited.more", "cited.once"), c("cited.once", "other"))

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

par(oma=c(1,0,0,0))

##########################################################################

genetypes=c("pc.cited.more", "pc.cited.once", "pc.other", "lnc.cited.more", "lnc.cited.once",  "lnc.other")

xpos.genetypes=c(1, 3.5, 6, 2,  4.5, 7)
names(xpos.genetypes)=genetypes

col.genetypes=rep(c("indianred", "steelblue"), each=3)
names(col.genetypes)=genetypes

col.genecat=c("indianred", "steelblue")
names(col.genecat)=c("pc", "lnc")

labels=c("a", "b")
names(labels)=c("tumor", "nontumor")

mtext.labels=c("tumor", "adjacent tissue")
names(mtext.labels)=c("tumor", "nontumor")

##########################################################################

## expression levels in meric, tumor samples

for(tissue in c("tumor", "nontumor")){

    par(mar=c(2.5, 3.75, 2.5, 1.75))

    xlim=c(0.5,8)
    ylim=c(0,20)

    tinyy=0.25
    tinyx=0.25

    plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=ylim)

    ## expression in tumor samples

    for(type in genetypes){
        this.genes=get(type)
        this.xpos=xpos.genetypes[type]
        this.col=col.genetypes[type]
        this.exp=get(paste("meantpm", tissue, sep="."))

       ## segments(this.xpos, ylim[1], this.xpos, ylim[2]*0.9, lty=3, col=this.col)

        vioplot(log2(this.exp[this.genes]+1), h=1, add=T, axes=F, at=this.xpos, border="black", pchMed=21, colMed="black", colMed2="white", cex=0.95, col=this.col)
    }

    ## signs for the p-values

    for(cittype in c("cited.more", "cited.once", "other")){
        this.pc=get(paste("pc", cittype,sep="."))
        this.lnc=get(paste("lnc", cittype,sep="."))
        this.exp=get(paste("meantpm", tissue, sep="."))
        this.pval=wilcox.test(this.exp[this.pc], this.exp[this.lnc])$p.value

        this.xpos=xpos.genetypes[c(paste("pc", cittype,sep="."),paste("lnc", cittype,sep="."))]
        this.yrange=range(c(log2(this.exp[this.pc]+1)))
        this.smally=diff(this.yrange)/20
        this.ypos=this.yrange[2]+this.smally

        segments(this.xpos[1], this.ypos, this.xpos[2], this.ypos, col="gray40")

        text(format(this.pval, digits=1), x=mean(this.xpos), y=this.ypos+tinyy, adj=c(0.5,0.01), col="black", cex=0.9, xpd=NA)
    }


    ## pairwise comparisons

    smally.genecat=c(diff(ylim)/200, diff(ylim)/7.8)
    names(smally.genecat)=c("pc", "lnc")

    smally.comp=c(0, diff(ylim)/15)

    for(i in 1:2){
        this.type1=comparisons[[i]][1]
        this.type2=comparisons[[i]][2]

        for(genetype in c("lnc")){
            this.genes1=get(paste(genetype, this.type1, sep="."))
            this.genes2=get(paste(genetype, this.type2, sep="."))

            this.pval=wilcox.test(this.exp[this.genes1], this.exp[this.genes2])$p.value

            this.xpos=xpos.genetypes[c(paste(genetype, this.type1, sep="."), paste(genetype, this.type2, sep="."))]
            this.ypos=ylim[2]-smally.genecat[genetype]-smally.comp[i]

            segments(this.xpos[1], this.ypos, this.xpos[2], this.ypos, col=col.genecat[genetype])

            text(format(this.pval, digits=1), x=mean(this.xpos), y=this.ypos+tinyy, adj=c(0.5,0.01), col=col.genecat[genetype], cex=0.9, xpd=NA)
        }
    }

    mtext(mtext.labels[tissue], side=3, line=0.5, at=mean(xpos.genetypes), cex=0.75)

    axis(side=2, mgp=c(3,0.65,0))
    mtext("mean expression level (log2 TPM)", side=2, line=2.5, cex=0.75)

    mtext("cited", side=1, at=mean(xpos.genetypes[c("pc.cited.more", "lnc.cited.more")]),line=0.5, cex=0.75)
    mtext(">1 articles", side=1, at=mean(xpos.genetypes[c("pc.cited.more", "lnc.cited.more")]),line=1.5, cex=0.75)

    mtext("cited", side=1, at=mean(xpos.genetypes[c("pc.cited.once", "lnc.cited.once")]),line=0.5, cex=0.75)
    mtext("1 article", side=1, at=mean(xpos.genetypes[c("pc.cited.once", "lnc.cited.once")]),line=1.5, cex=0.75)

    mtext("not cited", side=1, at=mean(xpos.genetypes[c("pc.other", "lnc.other")]),line=1, cex=0.75)

    axis(side=1, at=xpos.genetypes, labels=rep("", length(xpos.genetypes)), mgp=c(3,0.5,0))


    mtext(labels[tissue], font=2, line=0.95, at=-1)

}

##########################################################################

## antisense overlap and bidirectional promoters

labels=c("c", "d")
names(labels)=c("antisense.overlaps", "biprom")

mtext.labels=c("% with antisense overlap", "% with bidirectional promoters")
names(mtext.labels)=c("antisense.overlaps", "biprom")

smally.analysis=c(20, 35)
names(smally.analysis)=c("antisense.overlaps", "biprom")

for(analysis in c("antisense.overlaps", "biprom")){

    this.list.prop=get(paste("prop", analysis, sep="."))
    this.list.nb=get(paste("nb", analysis, sep="."))
    this.list.conf=get(paste(analysis, "conf", sep="."))

    par(mar=c(2.5, 3.75, 3, 1.75))

    ylim=c(0, 100)
    tinyy=2
    tinyx=0.25

    plot(1, type="n", xlab="", ylab="", axes=F, ylim=ylim, xlim=c(0.5,8))

    width=diff(xpos.genetypes)[1]/8

    for(type in genetypes){
        this.prop=100*this.list.prop[[type]]
        this.xpos=xpos.genetypes[type]
        this.col=col.genetypes[type]

        this.conf=this.list.conf[[type]]

        rect(this.xpos-width, 0, this.xpos+width, this.prop, col=this.col)
        segments(this.xpos, 100*this.conf[1], this.xpos, 100*this.conf[2], lwd=1.5)
    }

    mtext("cited", side=1, at=mean(xpos.genetypes[c("pc.cited.more", "lnc.cited.more")]),line=0.75, cex=0.75)
    mtext(">1", side=1, at=mean(xpos.genetypes[c("pc.cited.more", "lnc.cited.more")]),line=1.75, cex=0.75)

    mtext("cited", side=1, at=mean(xpos.genetypes[c("pc.cited.once", "lnc.cited.once")]),line=0.75, cex=0.75)
    mtext("1", side=1, at=mean(xpos.genetypes[c("pc.cited.once", "lnc.cited.once")]),line=1.75, cex=0.75)

    mtext("not cited", side=1, at=mean(xpos.genetypes[c("pc.other", "lnc.other")]),line=1.25, cex=0.75)

    axis(side=2, mgp=c(3,0.65,0))

    axis(side=1, at=xpos.genetypes, labels=rep("",length(xpos.genetypes)))

    ## abline(v=mean(xpos.genetypes[c("lnc.cited.more", "pc.cited.once")]), lty=3)
    ## abline(v=mean(xpos.genetypes[c("lnc.cited.once", "pc.other")]), lty=3)

    mtext(mtext.labels[analysis], side=2, line=2.5, cex=0.8)

    mtext(labels[analysis], font=2, line=0.95, at=-2)

    ## p-values, pc vs lnc

    for(cittype in c("cited.more", "cited.once", "other")){
        this.pc=get(paste("pc", cittype,sep="."))
        this.lnc=get(paste("lnc", cittype,sep="."))

        this.nb.pc.analysis=this.list.nb[[paste("pc", cittype, sep=".")]]
        this.nb.lnc.analysis=this.list.nb[[paste("lnc", cittype, sep=".")]]
        this.nb.pc.tot=length(this.pc)
        this.nb.lnc.tot=length(this.lnc)

        this.pval=prop.test(c(this.nb.pc.analysis, this.nb.lnc.analysis), c(this.nb.pc.tot, this.nb.lnc.tot))$p.value

        this.xpos=xpos.genetypes[c(paste("pc", cittype,sep="."),paste("lnc", cittype,sep="."))]

        this.ypos=ylim[2]-smally.analysis[analysis]

        segments(this.xpos[1]-tinyx, this.ypos, this.xpos[2]+tinyx, this.ypos, col="black")

        text(format(this.pval, digits=1), x=mean(this.xpos), y=this.ypos+tinyy, adj=c(0.5,0.01), col="black", cex=0.9)
    }

    ## p-values, pairwise comparisons

    smally.genecat=c(diff(ylim)/200, diff(ylim)/40.8)
    names(smally.genecat)=c("pc", "lnc")

    smally.comp=c(0, diff(ylim)/15)

    for(i in 1:2){
        this.type1=comparisons[[i]][1]
        this.type2=comparisons[[i]][2]

        for(genetype in c("lnc")){
            this.nb.type1.analysis=this.list.nb[[paste(genetype, this.type1, sep=".")]]
            this.nb.type2.analysis=this.list.nb[[paste(genetype, this.type2, sep=".")]]
            this.nb.type1.tot=length(get(paste(genetype, this.type1, sep=".")))
            this.nb.type2.tot=length(get(paste(genetype, this.type2, sep=".")))

            this.pval=prop.test(c(this.nb.type1.analysis, this.nb.type2.analysis), c(this.nb.type1.tot, this.nb.type2.tot))$p.value

            this.xpos=xpos.genetypes[c(paste(genetype, this.type1, sep="."), paste(genetype, this.type2, sep="."))]
            this.ypos=ylim[2]-smally.genecat[genetype]-smally.comp[i]

            segments(this.xpos[1], this.ypos, this.xpos[2], this.ypos, col=col.genecat[genetype])

            text(format(this.pval, digits=1), x=mean(this.xpos), y=this.ypos+tinyy, adj=c(0.5,0.01), col=col.genecat[genetype], cex=0.9, xpd=NA)
        }
    }
}

##########################################################################

## sequence conservation scores

par(mar=c(2.5, 3.25, 3, 0.75))

xlim=c(0.5,8)
ylim=c(0,1.2)

tinyy=0.05
tinyx=0.25

plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=ylim)

for(type in genetypes){
    this.genes=get(type)
    this.phast=phastcons[this.genes, "Score"]
    this.col=col.genetypes[type]
    this.xpos=xpos.genetypes[type]

    vioplot(this.phast, h=0.02, add=T, axes=F, at=this.xpos, border="black", pchMed=21, colMed="black", colMed2="white", cex=0.95, col=this.col)
}

 ## comparisons pc-lnc

for(cittype in c("cited.more", "cited.once", "other")){
    this.pc=get(paste("pc", cittype,sep="."))
    this.lnc=get(paste("lnc", cittype,sep="."))

    this.pval=wilcox.test(phastcons[this.pc, "Score"], phastcons[this.lnc, "Score"])$p.value

    this.xpos=xpos.genetypes[c(paste("pc", cittype,sep="."),paste("lnc", cittype,sep="."))]
    this.yrange=range(c(phastcons[this.pc, "Score"], phastcons[this.lnc, "Score"]), na.rm=T)
    this.smally=diff(this.yrange)/20
    this.ypos=this.yrange[2]+this.smally

    segments(this.xpos[1], this.ypos, this.xpos[2], this.ypos, col="gray40")

    text(format(this.pval, digits=1), x=mean(this.xpos), y=this.ypos+tinyy, adj=c(0.5,0.01), col="black", cex=0.9, xpd=NA)
}

## comparisons categories

smally.genecat=c(diff(ylim)/200, diff(ylim)/40.8)
names(smally.genecat)=c("pc", "lnc")

smally.comp=c(-diff(ylim)/10, 0)

for(i in 1:2){
    this.type1=comparisons[[i]][1]
    this.type2=comparisons[[i]][2]

    for(genetype in c("lnc")){
        this.genes1=get(paste(genetype, this.type1, sep="."))
        this.genes2=get(paste(genetype, this.type2, sep="."))

        this.pval=wilcox.test(phastcons[this.genes1, "Score"], phastcons[this.genes2, "Score"])$p.value

        this.xpos=xpos.genetypes[c(paste(genetype, this.type1, sep="."), paste(genetype, this.type2, sep="."))]
        this.ypos=ylim[2]-smally.genecat[genetype]-smally.comp[i]

        segments(this.xpos[1], this.ypos, this.xpos[2], this.ypos, col=col.genecat[genetype], xpd=NA)

        text(format(this.pval, digits=1), x=mean(this.xpos), y=this.ypos+tinyy, adj=c(0.5,0.01), col=col.genecat[genetype], cex=0.9, xpd=NA)
    }
}

axis(side=2, mgp=c(3,0.65,0), at=seq(from=0, to=1, by=0.2))


axis(side=1, at=xpos.genetypes, labels=rep("",length(xpos.genetypes)))
mtext("sequence conservation", side=2, line=2.5, cex=0.8, at=0.5)


mtext("cited", side=1, at=mean(xpos.genetypes[c("pc.cited.more", "lnc.cited.more")]),line=0.75, cex=0.75)
mtext(">1", side=1, at=mean(xpos.genetypes[c("pc.cited.more", "lnc.cited.more")]),line=1.75, cex=0.75)

mtext("cited", side=1, at=mean(xpos.genetypes[c("pc.cited.once", "lnc.cited.once")]),line=0.75, cex=0.75)
mtext("1", side=1, at=mean(xpos.genetypes[c("pc.cited.once", "lnc.cited.once")]),line=1.75, cex=0.75)

mtext("not cited", side=1, at=mean(xpos.genetypes[c("pc.other", "lnc.other")]),line=1.25, cex=0.75)

mtext("e", font=2, line=0.95, at=-1.35)

##########################################################################

legend("bottomleft", legend=c("protein-coding", "lncRNAs"), fill=c("indianred", "steelblue"), xpd=NA, inset=c(-2.25, -0.4), cex=1.1, box.col="white", bg="white", horiz=T)

##########################################################################

dev.off()

##########################################################################
