##########################################################################

if(!("pathFigures" %in% objects())){
    pathRData="../../data_for_publication/RData/"
    pathFigures="../../data_for_publication/figures/"

    maxFDR=0.05
    minLFC=0

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
    lnc.other=setdiff(lnc, lnc.cited.all)

    pc.cited.once=names(nb.citations.pc)[which(nb.citations.pc==1)]
    pc.cited.more=names(nb.citations.pc)[which(nb.citations.pc>1)]
    pc.cited.all=names(nb.citations.pc)
    pc.other=setdiff(pc, pc.cited.all)

    ## proportion differentially expressed

    nb.neighbor.significant=list() ## in what fraction of cases the gene has a close bidirectional promoter and the neighbor is significant

    for(type in c("tnt.meric", "grades", "tnt.tcga")){
        this.diffexp=get(paste("diffexp",type,sep="."))

        nb.neighbor.significant[[type]]=list()

        this.biprom=biprom1kb

        signif.genes=rownames(this.diffexp)[which(this.diffexp$padj < maxFDR)]

        this.biprom$NeighborSignificant=unlist(lapply(this.biprom$GenesCloseTSS, function(x) any(unlist(strsplit(x, split=","))%in%signif.genes)))

        for(genetype in c("lnc.cited.once", "lnc.cited.more", "lnc.other", "pc.cited.once", "pc.cited.more", "pc.other")){
            genes=get(genetype)

            nb.neighbor.significant[[type]][[genetype]]=length(which(genes%in%this.biprom$GeneID[which(this.biprom$NeighborSignificant)]))
        }
    }

    prepare=FALSE
}

##########################################################################

pdf(paste(pathFigures, "SupplementaryFigure7.pdf", sep=""), width=6.85, height=3)

##########################################################################

m=matrix(rep(NA, 1*15), nrow=1)

m[1,]=c(rep(1,5), rep(2,5), rep(3, 5))

layout(m)

##########################################################################

legends=c("tumor vs. adjacent tissue\nMERiC", "tumor grades\nMERiC", "tumor vs. adjacent tissue\nTCGA")
names(legends)=c("tnt.meric", "grades", "tnt.tcga")

genetypes=c("pc.cited.more", "pc.cited.once", "pc.other", "lnc.cited.more", "lnc.cited.once",  "lnc.other")

lty.genetypes=rep(1:3,2)
names(lty.genetypes)=genetypes

xpos.genetypes=c(1, 3.5, 6, 2,  4.5, 7)
names(xpos.genetypes)=genetypes

width=diff(xpos.genetypes)[1]/8

col.genetypes=rep(c("indianred", "steelblue"), each=3)
names(col.genetypes)=genetypes

## pairwise comparisons

comparisons=list(c("cited.more", "cited.once"), c("cited.once", "other"))


labels=letters[1:3]

##########################################################################

plotindex=1

for(type in c("tnt.meric", "grades", "tnt.tcga")){

    par(mar=c(3, 3.75, 3.5, 1))

    if(type=="tnt.meric"){
        ylim=c(0,50)
        tinyy=2
    }

    if(type=="grades"){
        ylim=c(0,30)
        tinyy=1.25
     }

    if(type=="tnt.tcga"){
        ylim=c(0,50)
        tiny=2
    }


    ## graphic pars

    diffy.comparisons=c(diff(ylim)/40, diff(ylim)/12)
    tinyx=0.25


    plot(1, type="n", xlab="", ylab="", axes=F, ylim=ylim, xlim=c(0.5,8))

    for(genetype in genetypes){
        genes=get(genetype)

        nb.signif.neighbors=nb.neighbor.significant[[type]][[genetype]]
        nb.tot=length(genes)

        this.prop=100*nb.signif.neighbors/nb.tot
        this.conf=100*prop.test(nb.signif.neighbors, nb.tot)$conf

        this.xpos=xpos.genetypes[genetype]
        this.col=col.genetypes[genetype]

        rect(this.xpos-width, 0, this.xpos+width, this.prop, col=this.col)
        segments(this.xpos, this.conf[1], this.xpos, this.conf[2], lwd=1.5)
    }


    ## p-values, comparison between lnc and pc, for each citation class

    for(citation in c("cited.more", "cited.once", "other")){
        nb.signif.pc=nb.neighbor.significant[[type]][[paste("pc", citation, sep=".")]]
        pc=get(paste("pc", citation, sep="."))
        nb.tot.pc=length(pc)

        nb.signif.lnc=nb.neighbor.significant[[type]][[paste("lnc", citation, sep=".")]]
        pc=get(paste("lnc", citation, sep="."))
        nb.tot.lnc=length(lnc)

        this.pval=prop.test(c(nb.signif.pc, nb.signif.lnc), c(nb.tot.pc, nb.tot.lnc))$p.value

        this.xpos=xpos.genetypes[c(paste("pc", citation,sep="."),paste("lnc", citation,sep="."))]
        this.ypos=ylim[2]*0.75

        segments(this.xpos[1]-tinyx, this.ypos, this.xpos[2]+tinyx, this.ypos, col="black")

        text(format(this.pval, digits=1), x=mean(this.xpos), y=this.ypos+tinyy, adj=c(0.5,0.01), col="black", cex=0.9)
    }


    ## p-values, comparison between lnc classes
    for(i in 1:2){
        this.type1=comparisons[[i]][1]
        this.type2=comparisons[[i]][2]

        for(genetype in c("lnc")){
            nb.signif.1=nb.neighbor.significant[[type]][[paste(genetype, this.type1, sep=".")]]
            genes1=get(paste(genetype, this.type1, sep="."))
            nb.tot.1=length(genes1)

            nb.signif.2=nb.neighbor.significant[[type]][[paste(genetype, this.type2, sep=".")]]
            genes2=get(paste(genetype, this.type2, sep="."))
            nb.tot.2=length(genes2)

            this.pval=prop.test(c(nb.signif.1, nb.signif.2), c(nb.tot.1, nb.tot.2))$p.value
            this.xpos=xpos.genetypes[c(paste(genetype, this.type1,sep="."),paste(genetype, this.type2,sep="."))]

            this.ypos=ylim[2]-diffy.comparisons[i]

            segments(this.xpos[1]-tinyx, this.ypos, this.xpos[2]+tinyx, this.ypos, col=col.genetypes[paste(genetype, this.type1, sep=".")], xpd=NA)

            text(format(this.pval, digits=1), x=mean(this.xpos), y=this.ypos+tinyy, adj=c(0.5,0.01), col=col.genetypes[paste(genetype, this.type1, sep=".")], cex=0.9, xpd=NA)
        }
    }


    mtext("cited", side=1, at=mean(xpos.genetypes[c("pc.cited.more", "lnc.cited.more")]),line=0.75, cex=0.75)
    mtext(">1 articles", side=1, at=mean(xpos.genetypes[c("pc.cited.more", "lnc.cited.more")]),line=1.75, cex=0.75)

    mtext("cited", side=1, at=mean(xpos.genetypes[c("pc.cited.once", "lnc.cited.once")]),line=0.75, cex=0.75)
    mtext("1 article", side=1, at=mean(xpos.genetypes[c("pc.cited.once", "lnc.cited.once")]),line=1.75, cex=0.75)

    mtext("not cited", side=1, at=mean(xpos.genetypes[c("pc.other", "lnc.other")]),line=1.25, cex=0.75)

    axis(side=2, mgp=c(3,0.65,0))

    axis(side=1, at=xpos.genetypes, labels=rep("",length(xpos.genetypes)))

    mtext("% with DE neighbors", side=2, line=2.5, cex=0.75)

    mtext(legends[type], line=1, at=4, side=3, cex=0.75)

    mtext(labels[plotindex], font=2, side=3, line=1.5, at=-1.65, cex=0.95)

    plotindex=plotindex+1
}

##########################################################################

dev.off()

##########################################################################
