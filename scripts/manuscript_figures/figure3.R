##########################################################################

if(!("pathFigures" %in% objects())){
    pathRData="../../data_for_publication/RData/"
    pathFigures="../../data_for_publication/main_figures/"

    maxFDR=0.01

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

    ## gene overlaps
    load(paste(pathRData, "data.gene.overlaps.RData",sep=""))

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
    nb.neighbor.significant=list() ## in what fraction of cases the gene has a close bidirectional promoter and the neighbor is signfiicant

    for(type in c("tissues", "grades", "tnt")){
        this.diffexp=get(paste("diffexp",type,sep="."))

        nb.significant[[type]]=list()
        nb.neighbor.significant[[type]]=list()

        this.biprom=biprom1kb

        signif.genes=rownames(this.diffexp)[which(this.diffexp$padj < maxFDR)]

        this.biprom$NeighborSignificant=unlist(lapply(this.biprom$GenesCloseTSS, function(x) any(unlist(strsplit(x, split=","))%in%signif.genes)))

        for(genetype in c("lnc.cited.once", "lnc.cited.more", "other.lnc", "pc.cited.once", "pc.cited.more", "other.pc")){
            genes=get(genetype)

            nb.significant[[type]][[genetype]]=length(which(this.diffexp$padj < maxFDR & rownames(this.diffexp)%in%genes))
            nb.neighbor.significant[[type]][[genetype]]=length(which(genes%in%this.biprom$GeneID[which(this.biprom$NeighborSignificant)]))
        }
    }

    prepare=FALSE
}

##########################################################################

pdf(paste(pathFigures, "Figure3.pdf", sep=""), width=5.85, height=7.1)

##########################################################################

m=matrix(rep(NA, 3*10), nrow=3)

m[1,]=c(rep(1,5), rep(2,5))
m[2,]=c(rep(3,5), rep(4,5))
m[3,]=c(rep(5,5), rep(6,5))

layout(m)

##########################################################################

legends=c("tumor vs. liver", "high grades vs. low grades", "tumor vs. adjacent tissue")
names(legends)=c("tissues", "grades", "tnt")

genetypes=c("pc.cited.more", "pc.cited.once", "other.pc", "lnc.cited.more", "lnc.cited.once",  "other.lnc")

xpos.genetypes=c(1, 3.5, 6, 2,  4.5, 7)
names(xpos.genetypes)=genetypes

width=diff(xpos.genetypes)[1]/8

col.genetypes=rep(c("indianred", "steelblue"), each=3)
names(col.genetypes)=genetypes

labels=letters[1:6]

##########################################################################

plotindex=1

for(type in c("tissues", "grades", "tnt")){

    ## proportion of significant genes

    par(mar=c(3, 3.75, 2.5, 1.75))

    plot(1, type="n", xlab="", ylab="", axes=F, ylim=c(0,60), xlim=c(0.5,8))

    for(genetype in genetypes){
        genes=get(genetype)

        nb.signif=nb.significant[[type]][[genetype]]
        nb.tot=length(genes)

        this.prop=100*nb.signif/nb.tot
        this.conf=100*prop.test(nb.signif, nb.tot)$conf

        this.xpos=xpos.genetypes[genetype]
        this.col=col.genetypes[genetype]

        rect(this.xpos-width, 0, this.xpos+width, this.prop, col=this.col)
        segments(this.xpos, this.conf[1], this.xpos, this.conf[2], lwd=1.5)
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

    mtext("% differentially expressed", side=2, line=2.5, cex=0.75)

    mtext(paste("differential expression,", legends[type]), line=1, at=9.5, side=3, cex=0.75)

    mtext(labels[plotindex], font=2, side=3, line=-0.5, at=-1.25, cex=0.85)
    plotindex=plotindex+1

    ##########################################################################

     ## proportion of significant neighbors

    par(mar=c(3, 3.75, 2.5, 1.75))

    plot(1, type="n", xlab="", ylab="", axes=F, ylim=c(0,30), xlim=c(0.5,8))

    for(genetype in genetypes){
        genes=get(genetype)

        nb.signif=nb.neighbor.significant[[type]][[genetype]]
        nb.tot=length(genes)

        this.prop=100*nb.signif/nb.tot
        this.conf=100*prop.test(nb.signif, nb.tot)$conf

        this.xpos=xpos.genetypes[genetype]
        this.col=col.genetypes[genetype]

        rect(this.xpos-width, 0, this.xpos+width, this.prop, col=this.col)
        segments(this.xpos, this.conf[1], this.xpos, this.conf[2], lwd=1.5)
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

    mtext("% DE neighbors", side=2, line=2.5, cex=0.75)

    mtext(labels[plotindex], font=2, side=3, line=-0.5, at=-1.25, cex=0.85)
    plotindex=plotindex+1

}

##########################################################################

dev.off()

##########################################################################
##########################################################################
