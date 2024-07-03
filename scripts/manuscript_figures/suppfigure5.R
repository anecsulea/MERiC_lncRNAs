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

    nb.significant=list()
    nb.upregulated=list()
    nb.downregulated=list()
    nb.total=list()

    for(type in c("tnt.tcga")){
        this.diffexp=get(paste("diffexp",type,sep="."))

        nb.significant[[type]]=list()
        nb.upregulated[[type]]=list()
        nb.downregulated[[type]]=list()

        signif.genes=rownames(this.diffexp)[which(this.diffexp$padj < maxFDR)]

        for(genetype in c("lnc.cited.once", "lnc.cited.more", "lnc.other", "pc.cited.once", "pc.cited.more", "pc.other")){
            genes=get(genetype)

            nb.total[[genetype]]=length(genes)

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

genetypes=c("pc.cited.more", "pc.cited.once", "pc.other", "lnc.cited.more", "lnc.cited.once",  "lnc.other")

lty.genetypes=rep(1:3,2)
names(lty.genetypes)=genetypes

xpos.genetypes=c(1, 3.5, 6, 2,  4.5, 7)
names(xpos.genetypes)=genetypes

width=diff(xpos.genetypes)[1]/8

col.genetypes=rep(c("indianred", "steelblue"), each=3)
names(col.genetypes)=genetypes

labels=letters[1:6]

## pairwise comparisons

comparisons=list(c("cited.more", "cited.once"), c("cited.once", "other"))

##########################################################################

## common graphic pars

par(mar=c(3, 3.75, 3.5, 1))

ylim=c(0,62)

diffy.comparisons=c(diff(ylim)/40, diff(ylim)/12)
tinyx=0.25
tinyy=2

##########################################################################

plotindex=1

for(type in c("tnt.tcga")){

    for(direction in c("up", "down")){

        this.nb.list=get(paste("nb.",direction,"regulated",sep=""))

        plot(1, type="n", xlab="", ylab="", axes=F, ylim=ylim, xlim=c(0.5,8))

        for(genetype in genetypes){
            genes=get(genetype)

            nb.signif=this.nb.list[[type]][[genetype]]
            nb.tot=length(genes)

            this.prop=100*nb.signif/nb.tot
            this.conf=100*prop.test(nb.signif, nb.tot)$conf

            this.xpos=xpos.genetypes[genetype]
            this.col=col.genetypes[genetype]

            rect(this.xpos-width, 0, this.xpos+width, this.prop, col=this.col)
            segments(this.xpos, this.conf[1], this.xpos, this.conf[2], lwd=1.5)
        }

        ## p-values, comparison between lnc and pc, for each citation class

        for(citation in c("cited.more", "cited.once", "other")){
            nb.signif.pc=this.nb.list[[type]][[paste("pc", citation, sep=".")]]
            nb.tot.pc=nb.total[[paste("pc", citation, sep=".")]]

            nb.signif.lnc=this.nb.list[[type]][[paste("lnc", citation, sep=".")]]
            nb.tot.lnc=nb.total[[paste("lnc", citation, sep=".")]]

            this.pval=prop.test(c(nb.signif.pc, nb.signif.lnc), c(nb.tot.pc, nb.tot.lnc))$p.value

            this.xpos=xpos.genetypes[c(paste("pc", citation,sep="."),paste("lnc", citation,sep="."))]
            this.ypos=ylim[2]*0.8

            segments(this.xpos[1]-tinyx, this.ypos, this.xpos[2]+tinyx, this.ypos, col="black")

            text(format(this.pval, digits=1), x=mean(this.xpos), y=this.ypos+tinyy, adj=c(0.5,0.01), col="black", cex=0.9)
        }

        ## p-values, comparison between lnc classes
        for(i in 1:2){
            this.type1=comparisons[[i]][1]
            this.type2=comparisons[[i]][2]

            for(genetype in c("lnc")){
                nb.signif.1=this.nb.list[[type]][[paste(genetype, this.type1, sep=".")]]
                nb.tot.1=nb.total[[paste(genetype, this.type1, sep=".")]]

                nb.signif.2=this.nb.list[[type]][[paste(genetype, this.type2, sep=".")]]
                nb.tot.2=nb.total[[paste(genetype, this.type2, sep=".")]]

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

        mtext("% up-regulated", side=2, line=2.5, cex=0.75)

        mtext(paste("up-regulated in", legends[type]), line=1, at=4, side=3, cex=0.75)

        mtext(labels[plotindex], font=2, side=3, line=1.5, at=-1.65, cex=0.95)
        plotindex=plotindex+1

    }
}

##########################################################################

##########################################################################

dev.off()

##########################################################################
##########################################################################
