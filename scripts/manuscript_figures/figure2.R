##########################################################################

if(!("pathFigures"%in%objects())){
    pathRData="../../data_for_publication/RData/"
    pathFigures="../../data_for_publication/main_figures/"

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

    load(paste(pathRData, "data.expression.levels.TCGA.RData",sep=""))
    tpm.tcga=tpm

    ## sample info
    load(paste(pathRData, "data.sample.info.RData",sep=""))
    liver.samples.meric=liver.samples
    tumor.samples.meric=tumor.samples

    load(paste(pathRData, "data.sample.info.TCGA.RData",sep=""))
    sampleinfo.tcga=sampleinfo
    tumor.samples.tcga=sampleinfo.tcga[which(sampleinfo$sample_type=="Tumor"),]
    nontumor.samples.tcga=sampleinfo.tcga[which(sampleinfo$sample_type=="Non-Tumor"),]

    ## sequence conservation
    load(paste(pathRData, "data.phastcons.30way.RData",sep=""))

    load=FALSE
}

##########################################################################

if(prepare){
    ## other lnc, not cited
    other.lnc=setdiff(lnc, all.cited.lnc)

    ## extract overlaps with pc genes
    overlaps.pc=gene.overlaps[which(gene.overlaps$NeighborID%in%pc),]

    ## sense and antisense overlaps
    sense.overlaps.pc=overlaps.pc[which(overlaps.pc$Type=="sense"),]
    antisense.overlaps.pc=overlaps.pc[which(overlaps.pc$Type=="antisense"),]

    prop.antisense.overlaps.cited.lnc=length(which(all.cited.lnc%in%antisense.overlaps.pc$GeneID))/length(all.cited.lnc)
    prop.antisense.overlaps.other.lnc=length(which(other.lnc%in%antisense.overlaps.pc$GeneID))/length(other.lnc)
    prop.antisense.overlaps.pc=length(which(pc%in%antisense.overlaps.pc$GeneID))/length(pc)

    antisense.conf.cited.lnc=prop.test(length(which(all.cited.lnc%in%antisense.overlaps.pc$GeneID)), length(all.cited.lnc))$conf
    antisense.conf.other.lnc=prop.test(length(which(other.lnc%in%antisense.overlaps.pc$GeneID)), length(other.lnc))$conf
    antisense.conf.pc=prop.test(length(which(pc%in%antisense.overlaps.pc$GeneID)), length(pc))$conf

    ## bidirectional promoters with pc genes

    biprom1pc=biprom1kb[which(biprom1kb$GenesCloseTSS%in%pc),]
    biprom5pc=biprom1kb[which(biprom5kb$GenesCloseTSS%in%pc),]

    prop.biprom.cited.lnc=length(which(all.cited.lnc%in%biprom1pc$GeneID))/length(all.cited.lnc)
    prop.biprom.other.lnc=length(which(other.lnc%in%biprom1pc$GeneID))/length(other.lnc)
    prop.biprom.pc=length(which(pc%in%biprom1pc$GeneID))/length(pc)

    biprom.conf.cited.lnc=prop.test(length(which(all.cited.lnc%in%biprom1pc$GeneID)), length(all.cited.lnc))$conf
    biprom.conf.other.lnc=prop.test(length(which(other.lnc%in%biprom1pc$GeneID)), length(other.lnc))$conf
    biprom.conf.pc=prop.test(length(which(pc%in%biprom1pc$GeneID)), length(pc))$conf

    ## average expression level across various samples
    meantpm.liver.meric=apply(tpm.meric[ ,liver.samples.meric$biopsyID],1, mean)
    meantpm.tumor.meric=apply(tpm.meric[ ,tumor.samples.meric$tumor_biopsyID],1, mean)
    meantpm.tumor.tcga=apply(tpm.tcga[, tumor.samples.tcga$id], 1, mean)
    meantpm.nontumor.tcga=apply(tpm.tcga[, nontumor.samples.tcga$id], 1, mean)

    ## phastcons score

    phast.cited.lnc=phastcons[all.cited.lnc, "Score"]
    phast.other.lnc=phastcons[other.lnc, "Score"]
    phast.pc=phastcons[pc, "Score"]

    prepare=FALSE
}

##########################################################################

pdf(paste(pathFigures, "Figure2.pdf", sep=""), width=6.85, height=5.5)

##########################################################################

## layout

m=matrix(rep(NA, 2*20), nrow=2)

m[1,]=c(rep(1, 10), rep(2,10))
m[2,]=c(rep(3, 5), rep(4, 5), rep(5,6), rep(6,4))

layout(m)

##########################################################################

## expression levels in meric

par(mar=c(2.5, 3.75, 2.5, 1.75))

xlim=c(0.5,7.5)
ylim=c(0,20)

plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=ylim)

## expression in tumor samples

vioplot(log2(meantpm.tumor.meric[pc]+1), h=1, add=T, axes=F, at=1, border="black", pchMed=21, colMed="black", colMed2="white", cex=0.95, col="indianred")
vioplot(log2(meantpm.tumor.meric[all.cited.lnc]+1), h=1, add=T, axes=F, at=2, border="black", pchMed=21, colMed="black", colMed2="white", cex=0.95, col="steelblue")
vioplot(log2(meantpm.tumor.meric[other.lnc]+1), h=1, add=T, axes=F, at=3, border="black", pchMed=21, colMed="black", colMed2="white", cex=0.95, col="lightblue")

mtext("tumor", side=3, line=-2, at=2, cex=0.8)


vioplot(log2(meantpm.liver.meric[pc]+1), h=1, add=T, axes=F, at=5, border="black", pchMed=21, colMed="black", colMed2="white", cex=0.95, col="indianred")
vioplot(log2(meantpm.liver.meric[all.cited.lnc]+1), h=1, add=T, axes=F, at=6, border="black", pchMed=21, colMed="black", colMed2="white", cex=0.95, col="steelblue")
vioplot(log2(meantpm.liver.meric[other.lnc]+1), h=1, add=T, axes=F, at=7, border="black", pchMed=21, colMed="black", colMed2="white", cex=0.95, col="lightblue")


mtext("liver", side=3, line=-2, at=6, cex=0.8)


axis(side=2, mgp=c(3,0.65,0))
mtext("mean expression level", side=2, line=2.5, cex=0.8)

axis(side=1, at=c(1:3, 5:7), labels=rep("", 6), mgp=c(3,0.5,0))

mtext("MERiC dataset", side=3, line=0.5, cex=0.8)


## legend("bottomleft", legend=c("protein-coding"), fill=c("indianred"), xpd=NA, inset=c(-0.1, -0.25), horiz=T, cex=1.2, bty="n")
## legend("bottomleft", legend=c( "HCC-associated lncRNAs", "other lncRNAs"), fill=c("steelblue", "lightblue"), xpd=NA, inset=c(0.38, -0.25), horiz=T, cex=1.2, bty="n")

mtext("a", font=2, line=0.95, at=-0.85)

##########################################################################

## expression levels in TCGA

par(mar=c(2.5, 3.75, 2.5, 1.75))

xlim=c(0.5,7.5)
ylim=c(0,20)

plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=ylim)

## expression in tumor samples

vioplot(log2(meantpm.tumor.tcga[pc]+1), h=1, add=T, axes=F, at=1, border="black", pchMed=21, colMed="black", colMed2="white", cex=0.95, col="indianred")
vioplot(log2(meantpm.tumor.tcga[all.cited.lnc]+1), h=1, add=T, axes=F, at=2, border="black", pchMed=21, colMed="black", colMed2="white", cex=0.95, col="steelblue")
vioplot(log2(meantpm.tumor.tcga[other.lnc]+1), h=1, add=T, axes=F, at=3, border="black", pchMed=21, colMed="black", colMed2="white", cex=0.95, col="lightblue")

mtext("tumor", side=3, line=-2, at=2, cex=0.8)

vioplot(log2(meantpm.nontumor.tcga[pc]+1), h=1, add=T, axes=F, at=5, border="black", pchMed=21, colMed="black", colMed2="white", cex=0.95, col="indianred")
vioplot(log2(meantpm.nontumor.tcga[all.cited.lnc]+1), h=1, add=T, axes=F, at=6, border="black", pchMed=21, colMed="black", colMed2="white", cex=0.95, col="steelblue")
vioplot(log2(meantpm.nontumor.tcga[other.lnc]+1), h=1, add=T, axes=F, at=7, border="black", pchMed=21, colMed="black", colMed2="white", cex=0.95, col="lightblue")


mtext("adjacent tissue", side=3, line=-2, at=6, cex=0.8)


axis(side=2, mgp=c(3,0.65,0))
mtext("mean expression level", side=2, line=2.5, cex=0.8)

axis(side=1, at=c(1:3, 5:7), labels=rep("", 6), mgp=c(3,0.5,0))

mtext("TCGA dataset", side=3, line=0.5, cex=0.8)

mtext("b", font=2, line=0.95, at=-0.85)

##########################################################################

## antisense overlap with pc genes

par(mar=c(2.5, 3.75, 3, 1.75))

values=100*c(prop.antisense.overlaps.pc, prop.antisense.overlaps.cited.lnc, prop.antisense.overlaps.other.lnc)

b=barplot(values, space=1, col=c("indianred", "steelblue", "lightblue"), ylim=c(0,50), axes=F, xlim=c(0.5,6.5))
segments(b[1], 100*antisense.conf.pc[1], b[1], 100*antisense.conf.pc[2], lwd=1.5)
segments(b[2], 100*antisense.conf.cited.lnc[1], b[2], 100*antisense.conf.cited.lnc[2], lwd=1.5)
segments(b[3], 100*antisense.conf.other.lnc[1], b[3], 100*antisense.conf.other.lnc[2], lwd=1.5)

axis(side=2, mgp=c(3,0.65,0))

axis(side=1, at=b, labels=rep("",3))

mtext("% with antisense overlap", side=2, line=2.5, cex=0.8)

mtext("c", font=2, line=0.95, at=-2.25)

##########################################################################

## antisense overlap with pc genes

par(mar=c(2.5, 3.25, 3, 2.25))

values=100*c(prop.biprom.pc, prop.biprom.cited.lnc, prop.biprom.other.lnc)

b=barplot(values, space=1, col=c("indianred", "steelblue", "lightblue"), ylim=c(0,40), axes=F, xlim=c(0.5,6.5))
segments(b[1], 100*biprom.conf.pc[1], b[1], 100*biprom.conf.pc[2], lwd=1.5)
segments(b[2], 100*biprom.conf.cited.lnc[1], b[2], 100*biprom.conf.cited.lnc[2], lwd=1.5)
segments(b[3], 100*biprom.conf.other.lnc[1], b[3], 100*biprom.conf.other.lnc[2], lwd=1.5)

axis(side=2, mgp=c(3,0.65,0))

axis(side=1, at=b, labels=rep("",3))

mtext("% with bidirectional promoters", side=2, line=2.5, cex=0.8)

mtext("d", font=2, line=0.95, at=-2.25)

##########################################################################

## sequence conservation scores

par(mar=c(2.5, 3.25, 3, 2.25))

xlim=c(0.5,3.5)
ylim=c(0,1)

plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=ylim)


vioplot(phast.pc, h=0.02, add=T, axes=F, at=1, border="black", pchMed=21, colMed="black", colMed2="white", cex=0.95, col="indianred")
vioplot(phast.cited.lnc, h=0.02, add=T, axes=F, at=2, border="black", pchMed=21, colMed="black", colMed2="white", cex=0.95, col="steelblue")
vioplot(phast.other.lnc, h=0.02, add=T, axes=F, at=3, border="black", pchMed=21, colMed="black", colMed2="white", cex=0.95, col="lightblue")

axis(side=2, mgp=c(3,0.65,0))

axis(side=1, at=1:3, labels=rep("",3))

mtext("sequence conservation", side=2, line=2.5, cex=0.8)

mtext("e", font=2, line=0.95, at=-0.6)

##########################################################################

## legend plot


par(mar=c(2.5, 0, 3, 1.75))


plot(1, type="n", xlab="", ylab="", axes=F)

legend("topleft", legend=c("protein-coding", "HCC-associated lncRNAs", "other lncRNAs"), fill=c("indianred","steelblue", "lightblue"), xpd=NA, inset=c(-0.4, -0.1), cex=1.1, bty="n")


##########################################################################

dev.off()

##########################################################################
