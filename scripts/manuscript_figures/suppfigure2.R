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


    load(paste(pathRData, "data.expression.levels.TCGA.RData",sep=""))
    tpm.tcga=tpm

    ## sample info

    load(paste(pathRData, "data.sample.info.TCGA.RData",sep=""))
    sampleinfo.tcga=sampleinfo
    tumor.samples.tcga=sampleinfo.tcga[which(sampleinfo$sample_type=="Tumor"),]
    nontumor.samples.tcga=sampleinfo.tcga[which(sampleinfo$sample_type=="Non-Tumor"),]

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

    ## average expression level across various samples

    meantpm.tumor.tcga=apply(tpm.tcga[, tumor.samples.tcga$id], 1, mean)
    meantpm.nontumor.tcga=apply(tpm.tcga[, nontumor.samples.tcga$id], 1, mean)

    prepare=FALSE
}

##########################################################################

pdf(paste(pathFigures, "SupplementaryFigure2.pdf", sep=""), width=5.85, height=3.25)

##########################################################################

## layout

m=matrix(rep(NA, 1*20), nrow=1)

m[1,]=c(rep(1, 10), rep(2,10))

layout(m)

par(oma=c(2,0,0,0))

##########################################################################

genetypes=c("pc.cited.more", "pc.cited.once", "pc.other", "lnc.cited.more", "lnc.cited.once",  "lnc.other")

xpos.genetypes=c(1, 3.5, 6, 2,  4.5, 7)
names(xpos.genetypes)=genetypes

col.genetypes=rep(c("indianred", "steelblue"), each=3)
names(col.genetypes)=genetypes

col.genecat=c("indianred", "steelblue")
names(col.genecat)=c("pc", "lnc")

tinyy=0.25
tinyx=0.25

labels=c("a", "b")
names(labels)=c("tumor", "nontumor")

titles=c("tumor", "adjacent tissue")
names(titles)=c("tumor", "nontumor")

comparisons=list(c("cited.more", "cited.once"), c("cited.once", "other"))

##########################################################################

for(tissue in c("tumor", "nontumor")){

    this.exp=get(paste("meantpm", tissue, "tcga", sep="."))

    par(mar=c(4, 3.75, 2.5, 1.75))

    xlim=c(0.5,8)
    ylim=c(0,20)

    plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=ylim)

    ## expression in tumor samples

    for(type in genetypes){
        this.genes=get(type)
        this.xpos=xpos.genetypes[type]
        this.col=col.genetypes[type]

        vioplot(log2(this.exp[this.genes]+1), h=1, add=T, axes=F, at=this.xpos, border="black", pchMed=21, colMed="black", colMed2="white", cex=0.95, col=this.col)
    }


    ###########################

    ## p-values, comparison between pc and lnc

    for(cittype in c("cited.more", "cited.once", "other")){
        this.pc=get(paste("pc", cittype,sep="."))
        this.lnc=get(paste("lnc", cittype,sep="."))

        this.pval=wilcox.test(this.exp[this.pc], this.exp[this.lnc])$p.value

        this.xpos=xpos.genetypes[c(paste("pc", cittype,sep="."),paste("lnc", cittype,sep="."))]
        this.yrange=range(c(log2(this.exp[this.pc]+1)))
        this.smally=diff(this.yrange)/20
        this.ypos=this.yrange[2]+this.smally

        segments(this.xpos[1], this.ypos, this.xpos[2], this.ypos, col="gray40")

        text(format(this.pval, digits=1), x=mean(this.xpos), y=this.ypos+tinyy, adj=c(0.5,0.01), col="black", cex=0.9, xpd=NA)
    }

    ## comparisons categories of genes

    smally.genecat=c(diff(ylim)/13, diff(ylim)/6.5)
    names(smally.genecat)=c("pc", "lnc")

    smally.comp=c(-diff(ylim)/10, 0)

    for(i in 1:2){
        this.type1=comparisons[[i]][1]
        this.type2=comparisons[[i]][2]

        for(genetype in c("lnc")){
            this.genes1=get(paste(genetype, this.type1, sep="."))
            this.genes2=get(paste(genetype, this.type2, sep="."))

            this.pval=wilcox.test(this.exp[this.genes1], this.exp[this.genes2])$p.value

            this.xpos=xpos.genetypes[c(paste(genetype, this.type1, sep="."), paste(genetype, this.type2, sep="."))]
            this.ypos=ylim[2]-smally.genecat[genetype]-smally.comp[i]

            segments(this.xpos[1], this.ypos, this.xpos[2], this.ypos, col=col.genecat[genetype], xpd=NA)

            text(format(this.pval, digits=1), x=mean(this.xpos), y=this.ypos+tinyy, adj=c(0.5,0.01), col=col.genecat[genetype], cex=0.9, xpd=NA)
        }
    }

    ###########################

    mtext(titles[tissue], side=3, line=0.25, at=mean(xpos.genetypes), cex=0.75)

    axis(side=2, mgp=c(3,0.65,0))
    mtext("mean expression level (log2 TPM)", side=2, line=2.5, cex=0.75)

    mtext("cited", side=1, at=mean(xpos.genetypes[c("pc.cited.more", "lnc.cited.more")]),line=0.5, cex=0.75)
    mtext(">1 articles", side=1, at=mean(xpos.genetypes[c("pc.cited.more", "lnc.cited.more")]),line=1.5, cex=0.75)

    mtext("cited", side=1, at=mean(xpos.genetypes[c("pc.cited.once", "lnc.cited.once")]),line=0.5, cex=0.75)
    mtext("1 article", side=1, at=mean(xpos.genetypes[c("pc.cited.once", "lnc.cited.once")]),line=1.5, cex=0.75)

    mtext("not cited", side=1, at=mean(xpos.genetypes[c("pc.other", "lnc.other")]),line=1, cex=0.75)

    axis(side=1, at=xpos.genetypes, labels=rep("", length(xpos.genetypes)), mgp=c(3,0.5,0))


    mtext(labels[tissue], font=2, line=0.95, at=-1.3)

   ##########################################################################

}

legend("bottomleft", legend=c("protein-coding", "lncRNAs"), fill=c("indianred", "steelblue"), xpd=NA, inset=c(-1.5, -0.35), cex=1.1, box.col="white", bg="white", horiz=T)

############################################################################

dev.off()

############################################################################
