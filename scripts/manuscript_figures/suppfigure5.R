##########################################################################

if(!("pathFigures" %in% objects())){
    pathRData="../../data_for_publication/RData/"
    pathFigures="../../data_for_publication/main_figures/"

    maxFDR=0.05
    minLFC=0

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

    ## differential expression
    load(paste(pathRData, "data.diffexp.RData",sep=""))

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

    ## proportion differentially expressed

    nb.significant=list()
    nb.upregulated=list()
    nb.downregulated=list()

    for(type in c("tnt.tcga")){
        this.diffexp=get(paste("diffexp",type,sep="."))

        nb.significant[[type]]=list()
        nb.upregulated[[type]]=list()
        nb.downregulated[[type]]=list()

        signif.genes=rownames(this.diffexp)[which(this.diffexp$padj < maxFDR)]

        for(genetype in c("lnc.cited.once", "lnc.cited.more", "other.lnc", "pc.cited.once", "pc.cited.more", "other.pc")){
            genes=get(genetype)

            nb.significant[[type]][[genetype]]=length(which(this.diffexp$padj < maxFDR & rownames(this.diffexp)%in%genes))

            nb.upregulated[[type]][[genetype]]=length(which(this.diffexp$padj < maxFDR & rownames(this.diffexp)%in%genes & this.diffexp$log2FoldChange>minLFC))
            nb.downregulated[[type]][[genetype]]=length(which(this.diffexp$padj < maxFDR & rownames(this.diffexp)%in%genes & this.diffexp$log2FoldChange< (-minLFC)))
        }
    }

    prepare=FALSE
}

##########################################################################

pdf(paste(pathFigures, "SupplementaryFigure5.pdf", sep=""), width=4.85, height=3)

##########################################################################

m=matrix(rep(NA, 1*10), nrow=1)

m[1,]=c(rep(1,5), rep(2,5))

layout(m)

##########################################################################

legends=c("tumor")
names(legends)=c("tnt.tcga")

genetypes=c("pc.cited.more", "pc.cited.once", "other.pc", "lnc.cited.more", "lnc.cited.once",  "other.lnc")

lty.genetypes=rep(1:3,2)
names(lty.genetypes)=genetypes

xpos.genetypes=c(1, 3.5, 6, 2,  4.5, 7)
names(xpos.genetypes)=genetypes

width=diff(xpos.genetypes)[1]/8

col.genetypes=rep(c("indianred", "steelblue"), each=3)
names(col.genetypes)=genetypes

labels=letters[1:6]

##########################################################################

plotindex=1

for(type in c("tnt.tcga")){

    ## proportion of up-regulated genes

    par(mar=c(3, 3.75, 3.5, 1))

    if(type=="tnt.tcga"){
        ylim=c(0,52)
    }

    plot(1, type="n", xlab="", ylab="", axes=F, ylim=ylim, xlim=c(0.5,8))

    for(genetype in genetypes){
        genes=get(genetype)

        nb.signif=nb.upregulated[[type]][[genetype]]
        nb.tot=length(genes)

        this.prop=100*nb.signif/nb.tot
        this.conf=100*prop.test(nb.signif, nb.tot)$conf

        this.xpos=xpos.genetypes[genetype]
        this.col=col.genetypes[genetype]

        rect(this.xpos-width, 0, this.xpos+width, this.prop, col=this.col)
        segments(this.xpos, this.conf[1], this.xpos, this.conf[2], lwd=1.5)
    }


    mtext("cited", side=1, at=mean(xpos.genetypes[c("pc.cited.more", "lnc.cited.more")]),line=0.75, cex=0.75)
    mtext(">1 articles", side=1, at=mean(xpos.genetypes[c("pc.cited.more", "lnc.cited.more")]),line=1.75, cex=0.75)

    mtext("cited", side=1, at=mean(xpos.genetypes[c("pc.cited.once", "lnc.cited.once")]),line=0.75, cex=0.75)
    mtext("1 article", side=1, at=mean(xpos.genetypes[c("pc.cited.once", "lnc.cited.once")]),line=1.75, cex=0.75)

    mtext("not cited", side=1, at=mean(xpos.genetypes[c("other.pc", "other.lnc")]),line=1.25, cex=0.75)


    axis(side=2, mgp=c(3,0.65,0))

    axis(side=1, at=xpos.genetypes, labels=rep("",length(xpos.genetypes)))

    abline(v=mean(xpos.genetypes[c("lnc.cited.more", "pc.cited.once")]), lty=3)
    abline(v=mean(xpos.genetypes[c("lnc.cited.once", "other.pc")]), lty=3)

    mtext("% up-regulated", side=2, line=2.5, cex=0.75)

    mtext(paste("up-regulated in", legends[type]), line=1, at=4, side=3, cex=0.75)

    mtext(labels[plotindex], font=2, side=3, line=1.5, at=-1.65, cex=0.95)
    plotindex=plotindex+1

##########################################################################

  ## proportion of down-regulated genes

    par(mar=c(3, 3.75, 3.5, 1))

    if(type=="tnt.meric"){
        ylim=c(0,60)
    }

     if(type=="grades"){
        ylim=c(0,30)
    }

    plot(1, type="n", xlab="", ylab="", axes=F, ylim=ylim, xlim=c(0.5,8))

    for(genetype in genetypes){
        genes=get(genetype)

        nb.signif=nb.downregulated[[type]][[genetype]]
        nb.tot=length(genes)

        this.prop=100*nb.signif/nb.tot
        this.conf=100*prop.test(nb.signif, nb.tot)$conf

        this.xpos=xpos.genetypes[genetype]
        this.col=col.genetypes[genetype]

        rect(this.xpos-width, 0, this.xpos+width, this.prop, col=this.col)
        segments(this.xpos, this.conf[1], this.xpos, this.conf[2], lwd=1.5)
    }


    mtext("cited", side=1, at=mean(xpos.genetypes[c("pc.cited.more", "lnc.cited.more")]),line=0.75, cex=0.75)
    mtext(">1 articles", side=1, at=mean(xpos.genetypes[c("pc.cited.more", "lnc.cited.more")]),line=1.75, cex=0.75)

    mtext("cited", side=1, at=mean(xpos.genetypes[c("pc.cited.once", "lnc.cited.once")]),line=0.75, cex=0.75)
    mtext("1 article", side=1, at=mean(xpos.genetypes[c("pc.cited.once", "lnc.cited.once")]),line=1.75, cex=0.75)

    mtext("not cited", side=1, at=mean(xpos.genetypes[c("other.pc", "other.lnc")]),line=1.25, cex=0.75)


    axis(side=2, mgp=c(3,0.65,0))

    axis(side=1, at=xpos.genetypes, labels=rep("",length(xpos.genetypes)))

    abline(v=mean(xpos.genetypes[c("lnc.cited.more", "pc.cited.once")]), lty=3)
    abline(v=mean(xpos.genetypes[c("lnc.cited.once", "other.pc")]), lty=3)

    mtext("% down-regulated", side=2, line=2.5, cex=0.75)

    mtext(paste("down-regulated in", legends[type]), line=1, at=4, side=3, cex=0.75)

    mtext(labels[plotindex], font=2, side=3, line=1.5, at=-1.65, cex=0.95)
    plotindex=plotindex+1

}

##########################################################################

dev.off()

##########################################################################
##########################################################################
