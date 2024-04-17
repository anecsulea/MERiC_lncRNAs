#############################################################################

if(!("pathFigures"%in%objects())){
    pathRData="../../data_for_publication/RData/"
    pathFigures="../../data_for_publication/main_figures/"

    load=TRUE
    prepare=TRUE
}

#############################################################################

if(load){
    ## cited lncRNAs
    load(paste(pathRData, "data.PubMed.analysis.RData",sep=""))

    ## gene info
    load(paste(pathRData, "data.gene.info.RData", sep=""))

    ## differential expression
    load(paste(pathRData, "data.diffexp.RData",sep=""))

    maxFDR=0.05
    minFC=1.5

    load=FALSE
}

#############################################################################

if(prepare){

    lnc.cited.once=names(nb.citations.lnc)[which(nb.citations.lnc==1)]
    lnc.cited.more=names(nb.citations.lnc)[which(nb.citations.lnc>1)]
    lnc.cited.all=names(nb.citations.lnc)
    other.lnc=setdiff(lnc, lnc.cited.all)

    pc.cited.once=names(nb.citations.pc)[which(nb.citations.pc==1)]
    pc.cited.more=names(nb.citations.pc)[which(nb.citations.pc>1)]
    pc.cited.all=names(nb.citations.pc)
    other.pc=setdiff(pc, pc.cited.all)

    prepare=FALSE
}

#############################################################################

pdf(file=paste(pathFigures, "SupplementaryFigure7.pdf", sep=""), width=6.85, height=9.5)

m=matrix(rep(NA, 3*10), nrow=3)

m[1,]=c(rep(1,5), rep(2,5))
m[2,]=c(rep(3,5), rep(4,5))
m[3,]=c(rep(5,5), rep(6,5))

layout(m)

#############################################################################

genetypes=c("pc.cited.more", "lnc.cited.more", "pc.cited.once", "lnc.cited.once", "other.pc", "other.lnc")

titles=c("protein-coding, cited >1 articles", "lncRNAs, cited >1 articles",  "protein-coding, cited 1 article", "lncRNAs, cited 1 article", "protein-coding, not cited", "lncRNAs, not cited")
names(titles)=genetypes

labels=letters[1:6]
names(labels)=genetypes

for(genetype in genetypes){
    this.genes=get(genetype)

    lfc.meric=diffexp.tnt.meric[this.genes, "log2FoldChange"]
    lfc.tcga=diffexp.tnt.tcga[this.genes, "log2FoldChange"]

    padj.meric=diffexp.tnt.meric[this.genes, "padj"]
    padj.tcga=diffexp.tnt.tcga[this.genes, "padj"]

    lim=range(c(lfc.meric, lfc.tcga), na.rm=T)

    par(mar=c(4.1, 4.1, 2.1, 1.1))
    plot(1, type="n", xlab="", ylab="", axes=F, xlim=lim, ylim=lim)

    consistent=which((lfc.meric*lfc.tcga)>0 & padj.meric<maxFDR & padj.tcga<maxFDR)
    inconsistent=which((lfc.meric*lfc.tcga)<0 & padj.meric<maxFDR & padj.tcga<maxFDR)
    other=setdiff(1:length(lfc.meric), c(consistent, inconsistent))

    points(lfc.meric[other], lfc.tcga[other], pch=21, bg="gray40")
    points(lfc.meric[consistent], lfc.tcga[consistent], pch=21, bg="seagreen")
    points(lfc.meric[inconsistent], lfc.tcga[inconsistent], pch=21, bg="darkred")

    axis(side=1, mgp=c(3, 0.5, 0))
    axis(side=2, mgp=c(3, 0.65, 0))
    box()

    mtext("log2 fold change, MERiC", side=1, line=2, cex=0.85)
    mtext("log2 fold change, TCGA", side=2, line=2, cex=0.85)

    abline(h=0, lty=3, col="gray40")
    abline(v=0, lty=3, col="gray40")
    abline(0, 1, lty=3, col="gray40")

    rho=round(cor(lfc.meric, lfc.tcga, use="complete.obs", method="spearman"), digits=2)

    text(paste("rho =",rho), x=lim[1], y=lim[2], adj=c(0.1, 0.9), cex=1.1)

    mtext(titles[genetype], side=3, cex=0.85)
    mtext(labels[genetype], side=3, cex=1.1, font=2, at=lim[1]-diff(lim)/5.8, line=0.5)

    pc.consistent=round(100*length(consistent)/length(this.genes), digits=0)
    pc.inconsistent=round(100*length(inconsistent)/length(this.genes), digits=0)
    pc.other=round(100*length(other)/length(this.genes), digits=0)

    legend("bottomright", legend=c(paste("consistent DE (", pc.consistent,"%)",sep=""),paste("contradictory DE (", pc.inconsistent,"%)",sep=""), paste("other (", pc.other,"%)",sep="")), pch=21, inset=0.01, pt.bg=c("seagreen", "darkred", "gray40"))

}

#############################################################################

dev.off()

#############################################################################
