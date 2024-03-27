##########################################################################

if(!("pathFigures" %in% objects())){
    pathRData="../../data_for_publication/RData/"
    pathFigures="../../data_for_publication/main_figures/"

    maxFDR=0.05

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
    ## other lnc, not cited
    other.lnc=setdiff(lnc, all.cited.lnc)

    ## significant genes for volcano plot

    diffexp.tissues.significant=diffexp.tissues[which(diffexp.tissues$padj<0.05),]
    diffexp.grades.significant=diffexp.grades[which(diffexp.grades$padj<0.05),]
    diffexp.tnt.significant=diffexp.tnt[which(diffexp.tnt$padj<0.05),]

    ## proportion up-regulated and down-regulated

    prop=list()

    conf.up=list()
    conf.down=list()

    for(type in c("tissues", "grades", "tnt")){
        this.diffexp=get(paste("diffexp",type,sep="."))

        prop.up.cited.lnc=length(which(this.diffexp[all.cited.lnc, "padj"]<maxFDR & this.diffexp[all.cited.lnc, "log2FoldChange"]>0))/length(all.cited.lnc)
        prop.up.other.lnc=length(which(this.diffexp[other.lnc, "padj"]<maxFDR & this.diffexp[other.lnc, "log2FoldChange"]>0))/length(other.lnc)
        prop.up.pc=length(which(this.diffexp[pc, "padj"]<maxFDR & this.diffexp[pc, "log2FoldChange"]>0))/length(pc)

        prop.down.cited.lnc=length(which(this.diffexp[all.cited.lnc, "padj"]<maxFDR & this.diffexp[all.cited.lnc, "log2FoldChange"]<0))/length(all.cited.lnc)
        prop.down.other.lnc=length(which(this.diffexp[other.lnc, "padj"]<maxFDR & this.diffexp[other.lnc, "log2FoldChange"]<0))/length(other.lnc)
        prop.down.pc=length(which(this.diffexp[pc, "padj"]<maxFDR & this.diffexp[pc, "log2FoldChange"]<0))/length(pc)

        prop[[type]]=c(prop.up.pc, prop.up.cited.lnc, prop.up.other.lnc, prop.down.pc, prop.down.cited.lnc, prop.down.other.lnc)

        conf.up[[type]]=list()
        conf.up[[type]][["pc"]]=prop.test(length(which(this.diffexp[pc, "padj"]<maxFDR & this.diffexp[pc, "log2FoldChange"]>0)), length(pc))$conf
        conf.up[[type]][["cited.lnc"]]=prop.test(length(which(this.diffexp[all.cited.lnc, "padj"]<maxFDR & this.diffexp[all.cited.lnc, "log2FoldChange"]>0)), length(all.cited.lnc))$conf
        conf.up[[type]][["other.lnc"]]=prop.test(length(which(this.diffexp[other.lnc, "padj"]<maxFDR & this.diffexp[other.lnc, "log2FoldChange"]>0)), length(other.lnc))$conf

        conf.down[[type]]=list()
        conf.down[[type]][["pc"]]=prop.test(length(which(this.diffexp[pc, "padj"]<maxFDR & this.diffexp[pc, "log2FoldChange"]<0)), length(pc))$conf
        conf.down[[type]][["cited.lnc"]]=prop.test(length(which(this.diffexp[all.cited.lnc, "padj"]<maxFDR & this.diffexp[all.cited.lnc, "log2FoldChange"]<0)), length(all.cited.lnc))$conf
        conf.down[[type]][["other.lnc"]]=prop.test(length(which(this.diffexp[other.lnc, "padj"]<maxFDR & this.diffexp[other.lnc, "log2FoldChange"]<0)), length(other.lnc))$conf
    }

    prepare=FALSE
}

##########################################################################

pdf(paste(pathFigures, "Figure3.pdf", sep=""), width=5, height=7.5)

##########################################################################

m=matrix(rep(NA, 3*9), nrow=3)

m[1,]=c(rep(1,4), rep(2,5))
m[2,]=c(rep(3,4), rep(4,5))
m[3,]=c(rep(5,4), rep(6,5))

layout(m)

##########################################################################


labels=letters[1:9]
i=1

###########################

legends=c("MERiC\ntumor vs. liver", "MERiC\nhigh grades vs. low grades", "TCGA\ntumor vs. adjacent tissue")
names(legends)=c("tissues", "grades", "tnt")

for(type in c("tissues", "grades", "tnt")){
    this.diffexp.signif=get(paste("diffexp", type, "significant",sep="."))

    ###########################

    ## violin plot of fold changes for significant genes

    par(mar=c(4.5, 3.5, 2.5, 2.1))

    fmax=max(abs(this.diffexp.signif$log2FoldChange))
    ylim=c(-fmax, fmax)

    plot(1, type="n", xlab="", ylab="", axes=F, xlim=c(0.5,3.5), ylim=ylim)

    vioplot(this.diffexp.signif[intersect(rownames(this.diffexp.signif), pc), "log2FoldChange"], h=0.25, add=T, axes=F, at=1, border="black", pchMed=21, colMed="black", colMed2="white", cex=0.95, col="indianred")
    vioplot(this.diffexp.signif[intersect(rownames(this.diffexp.signif), all.cited.lnc), "log2FoldChange"], h=0.25, add=T, axes=F, at=2, border="black", pchMed=21, colMed="black", colMed2="white", cex=0.95, col="steelblue")
    vioplot(this.diffexp.signif[intersect(rownames(this.diffexp.signif), other.lnc), "log2FoldChange"], h=0.25, add=T, axes=F, at=3, border="black", pchMed=21, colMed="black", colMed2="white", cex=0.95, col="lightblue")

    axis(side=2, mgp=c(3,0.6,0))
    mtext("log2 fold change", side=2, line=2, cex=0.75)

   ## axis(side=1, mgp=c(3,0.5,0), at=1:3, labels=rep("",3))

    mtext(labels[i], side=3, line=1.25, font=2, at=-0.35)
    i=i+1

###########################

    par(mar=c(3.5, 2.5, 2.5, 3.1))

    ## barplot for % up-regulated down-regulated

    this.prop=100*prop[[type]]

    ylim=c(0, max(this.prop)+10)

    b=barplot(this.prop, col=rep(c("indianred", "steelblue", "lightblue"), 2), space=c(1, rep(0.75,2), 1.5, rep(0.75,2)), xlim=c(1,11.5), ylim=ylim)

    this.conf.up=conf.up[[type]]
    segments(b[1], 100*this.conf.up[["pc"]][1], b[1], 100*this.conf.up[["pc"]][2], lwd=1.5)
    segments(b[2], 100*this.conf.up[["cited.lnc"]][1], b[2], 100*this.conf.up[["cited.lnc"]][2], lwd=1.5)
    segments(b[3], 100*this.conf.up[["other.lnc"]][1], b[3], 100*this.conf.up[["other.lnc"]][2], lwd=1.5)

    mtext("up-regulated", at=b[2], side=1, line=0.25, cex=0.7)

    abline(v=mean(b[3:4]), lty=3, col="gray40")

    this.conf.down=conf.down[[type]]
    segments(b[4], 100*this.conf.down[["pc"]][1], b[4], 100*this.conf.down[["pc"]][2], lwd=1.5)
    segments(b[5], 100*this.conf.down[["cited.lnc"]][1], b[5], 100*this.conf.down[["cited.lnc"]][2], lwd=1.5)
    segments(b[6], 100*this.conf.down[["other.lnc"]][1], b[6], 100*this.conf.down[["other.lnc"]][2], lwd=1.5)

    mtext("down-regulated", at=b[5], side=1, line=0.25, cex=0.7)

    mtext("% of genes", side=2, line=2.5, cex=0.75)

    mtext(labels[i], side=3, line=1.25, font=2, at=-1.55)
    i=i+1

    if(type=="tnt"){
        legend("topright", legend=c("protein-coding", "HCC-associated lncRNAs", "other lncRNAs"), fill=c("indianred","steelblue", "lightblue"), xpd=NA, inset=c(-0.1, -0.3), cex=1.1, bty="n")
    }

    mtext(legends[type], side=4, cex=0.75, line=1.5)
}

##########################################################################

dev.off()

##########################################################################
