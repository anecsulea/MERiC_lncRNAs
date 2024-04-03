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
    lnc.cited.once=names(nb.citations.lnc)[which(nb.citations.lnc==1)]
    lnc.cited.more=names(nb.citations.lnc)[which(nb.citations.lnc>1)]
    lnc.cited.all=names(nb.citations.lnc)

    ## other lnc, not cited
    other.lnc=setdiff(lnc, lnc.cited.all)

    ## extract overlaps with pc genes
    overlaps.pc=gene.overlaps[which(gene.overlaps$NeighborID%in%pc),]

    ## sense and antisense overlaps
    sense.overlaps.pc=overlaps.pc[which(overlaps.pc$Type=="sense"),]
    antisense.overlaps.pc=overlaps.pc[which(overlaps.pc$Type=="antisense"),]

    prop.antisense.overlaps.cited.once.lnc=length(which(lnc.cited.once%in%antisense.overlaps.pc$GeneID))/length(lnc.cited.once)
    prop.antisense.overlaps.cited.more.lnc=length(which(lnc.cited.more%in%antisense.overlaps.pc$GeneID))/length(lnc.cited.more)
    prop.antisense.overlaps.other.lnc=length(which(other.lnc%in%antisense.overlaps.pc$GeneID))/length(other.lnc)
    prop.antisense.overlaps.pc=length(which(pc%in%antisense.overlaps.pc$GeneID))/length(pc)

    antisense.conf.cited.once.lnc=prop.test(length(which(lnc.cited.once%in%antisense.overlaps.pc$GeneID)), length(lnc.cited.once))$conf
    antisense.conf.cited.more.lnc=prop.test(length(which(lnc.cited.more%in%antisense.overlaps.pc$GeneID)), length(lnc.cited.more))$conf
    antisense.conf.other.lnc=prop.test(length(which(other.lnc%in%antisense.overlaps.pc$GeneID)), length(other.lnc))$conf
    antisense.conf.pc=prop.test(length(which(pc%in%antisense.overlaps.pc$GeneID)), length(pc))$conf

    ## bidirectional promoters with pc genes

    biprom1pc=biprom1kb[which(biprom1kb$GenesCloseTSS%in%pc),]
    biprom5pc=biprom1kb[which(biprom5kb$GenesCloseTSS%in%pc),]

    prop.biprom.cited.once.lnc=length(which(lnc.cited.once%in%biprom1pc$GeneID))/length(lnc.cited.once)
    prop.biprom.cited.more.lnc=length(which(lnc.cited.more%in%biprom1pc$GeneID))/length(lnc.cited.more)
    prop.biprom.other.lnc=length(which(other.lnc%in%biprom1pc$GeneID))/length(other.lnc)
    prop.biprom.pc=length(which(pc%in%biprom1pc$GeneID))/length(pc)

    biprom.conf.cited.once.lnc=prop.test(length(which(lnc.cited.once%in%biprom1pc$GeneID)), length(lnc.cited.once))$conf
    biprom.conf.cited.more.lnc=prop.test(length(which(lnc.cited.more%in%biprom1pc$GeneID)), length(lnc.cited.more))$conf
    biprom.conf.other.lnc=prop.test(length(which(other.lnc%in%biprom1pc$GeneID)), length(other.lnc))$conf
    biprom.conf.pc=prop.test(length(which(pc%in%biprom1pc$GeneID)), length(pc))$conf

    ## average expression level across various samples
    meantpm.liver.meric=apply(tpm.meric[ ,liver.samples.meric$biopsyID],1, mean)
    meantpm.tumor.meric=apply(tpm.meric[ ,tumor.samples.meric$tumor_biopsyID],1, mean)
    meantpm.tumor.tcga=apply(tpm.tcga[, tumor.samples.tcga$id], 1, mean)
    meantpm.nontumor.tcga=apply(tpm.tcga[, nontumor.samples.tcga$id], 1, mean)

    ## phastcons score

    phast.cited.once.lnc=phastcons[lnc.cited.once, "Score"]
    phast.cited.more.lnc=phastcons[lnc.cited.more, "Score"]
    phast.other.lnc=phastcons[other.lnc, "Score"]
    phast.pc=phastcons[pc, "Score"]

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

## expression levels in meric

par(mar=c(2.5, 3.75, 2.5, 1.75))

xlim=c(0.5,9.5)
ylim=c(0,20)

plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=ylim)

## expression in tumor samples

vioplot(log2(meantpm.tumor.meric[pc]+1), h=1, add=T, axes=F, at=1, border="black", pchMed=21, colMed="black", colMed2="white", cex=0.95, col="indianred")
vioplot(log2(meantpm.tumor.meric[lnc.cited.more]+1), h=1, add=T, axes=F, at=2, border="black", pchMed=21, colMed="black", colMed2="white", cex=0.95, col="dodgerblue4")
vioplot(log2(meantpm.tumor.meric[lnc.cited.once]+1), h=1, add=T, axes=F, at=3, border="black", pchMed=21, colMed="black", colMed2="white", cex=0.95, col="steelblue")
vioplot(log2(meantpm.tumor.meric[other.lnc]+1), h=1, add=T, axes=F, at=4, border="black", pchMed=21, colMed="black", colMed2="white", cex=0.95, col="lightblue")

mtext("tumor", side=3, line=-2, at=2.5, cex=0.8)

vioplot(log2(meantpm.liver.meric[pc]+1), h=1, add=T, axes=F, at=6, border="black", pchMed=21, colMed="black", colMed2="white", cex=0.95, col="indianred")
vioplot(log2(meantpm.liver.meric[lnc.cited.more]+1), h=1, add=T, axes=F, at=7, border="black", pchMed=21, colMed="black", colMed2="white", cex=0.95, col="dodgerblue4")
vioplot(log2(meantpm.liver.meric[lnc.cited.once]+1), h=1, add=T, axes=F, at=8, border="black", pchMed=21, colMed="black", colMed2="white", cex=0.95, col="steelblue")
vioplot(log2(meantpm.liver.meric[other.lnc]+1), h=1, add=T, axes=F, at=9, border="black", pchMed=21, colMed="black", colMed2="white", cex=0.95, col="lightblue")

mtext("liver", side=3, line=-2, at=7.5, cex=0.8)

axis(side=2, mgp=c(3,0.65,0))
mtext("mean expression level", side=2, line=2.5, cex=0.8)

axis(side=1, at=c(1:4, 6:9), labels=rep("", 8), mgp=c(3,0.5,0))

mtext("MERiC dataset", side=3, line=0.5, cex=0.8)

mtext("a", font=2, line=0.95, at=-1.3)

##########################################################################

## expression levels in TCGA

par(mar=c(2.5, 3.75, 2.5, 1.75))

xlim=c(0.5,9.5)
ylim=c(0,20)

plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=ylim)

## expression in tumor samples

vioplot(log2(meantpm.tumor.tcga[pc]+1), h=1, add=T, axes=F, at=1, border="black", pchMed=21, colMed="black", colMed2="white", cex=0.95, col="indianred")
vioplot(log2(meantpm.tumor.tcga[lnc.cited.more]+1), h=1, add=T, axes=F, at=2, border="black", pchMed=21, colMed="black", colMed2="white", cex=0.95, col="dodgerblue4")
vioplot(log2(meantpm.tumor.tcga[lnc.cited.once]+1), h=1, add=T, axes=F, at=3, border="black", pchMed=21, colMed="black", colMed2="white", cex=0.95, col="steelblue")
vioplot(log2(meantpm.tumor.tcga[other.lnc]+1), h=1, add=T, axes=F, at=4, border="black", pchMed=21, colMed="black", colMed2="white", cex=0.95, col="lightblue")

mtext("tumor", side=3, line=-2, at=2.5, cex=0.8)

vioplot(log2(meantpm.nontumor.tcga[pc]+1), h=1, add=T, axes=F, at=6, border="black", pchMed=21, colMed="black", colMed2="white", cex=0.95, col="indianred")
vioplot(log2(meantpm.nontumor.tcga[lnc.cited.more]+1), h=1, add=T, axes=F, at=7, border="black", pchMed=21, colMed="black", colMed2="white", cex=0.95, col="dodgerblue4")
vioplot(log2(meantpm.nontumor.tcga[lnc.cited.once]+1), h=1, add=T, axes=F, at=8, border="black", pchMed=21, colMed="black", colMed2="white", cex=0.95, col="steelblue")
vioplot(log2(meantpm.nontumor.tcga[other.lnc]+1), h=1, add=T, axes=F, at=9, border="black", pchMed=21, colMed="black", colMed2="white", cex=0.95, col="lightblue")

mtext("adjacent tissue", side=3, line=-2, at=7.5, cex=0.8)

axis(side=2, mgp=c(3,0.65,0))
mtext("mean expression level", side=2, line=2.5, cex=0.8)

axis(side=1, at=c(1:4, 6:9), labels=rep("", 8), mgp=c(3,0.5,0))

mtext("TCGA dataset", side=3, line=0.5, cex=0.8)

mtext("b", font=2, line=0.95, at=-1.3)

##########################################################################

## antisense overlap with pc genes

par(mar=c(2.5, 3.75, 3, 1.75))

values=100*c(prop.antisense.overlaps.pc, prop.antisense.overlaps.cited.more.lnc, prop.antisense.overlaps.cited.once.lnc, prop.antisense.overlaps.other.lnc)

b=barplot(values, space=1, col=c("indianred", "dodgerblue4", "steelblue", "lightblue"), ylim=c(0,60), axes=F, xlim=c(0.5,8))
segments(b[1], 100*antisense.conf.pc[1], b[1], 100*antisense.conf.pc[2], lwd=1.5)
segments(b[2], 100*antisense.conf.cited.more.lnc[1], b[2], 100*antisense.conf.cited.more.lnc[2], lwd=1.5)
segments(b[3], 100*antisense.conf.cited.once.lnc[1], b[3], 100*antisense.conf.cited.once.lnc[2], lwd=1.5)
segments(b[4], 100*antisense.conf.other.lnc[1], b[4], 100*antisense.conf.other.lnc[2], lwd=1.5)

axis(side=2, mgp=c(3,0.65,0))

axis(side=1, at=b, labels=rep("",4))

mtext("% with antisense overlap", side=2, line=2.5, cex=0.8)

mtext("c", font=2, line=0.95, at=-2)

##########################################################################

## antisense overlap with pc genes

par(mar=c(2.5, 3.25, 3, 2.25))

values=100*c(prop.biprom.pc, prop.biprom.cited.more.lnc,prop.biprom.cited.once.lnc, prop.biprom.other.lnc)

b=barplot(values, space=1, col=c("indianred", "dodgerblue4", "steelblue", "lightblue"), ylim=c(0,40), axes=F, xlim=c(0.5, 8))
segments(b[1], 100*biprom.conf.pc[1], b[1], 100*biprom.conf.pc[2], lwd=1.5)
segments(b[2], 100*biprom.conf.cited.more.lnc[1], b[2], 100*biprom.conf.cited.more.lnc[2], lwd=1.5)
segments(b[3], 100*biprom.conf.cited.once.lnc[1], b[3], 100*biprom.conf.cited.once.lnc[2], lwd=1.5)
segments(b[4], 100*biprom.conf.other.lnc[1], b[4], 100*biprom.conf.other.lnc[2], lwd=1.5)

axis(side=2, mgp=c(3,0.65,0))

axis(side=1, at=b, labels=rep("",4))

mtext("% with bidirectional promoters", side=2, line=2.5, cex=0.8)

mtext("d", font=2, line=0.95, at=-2)

##########################################################################

## sequence conservation scores

par(mar=c(2.5, 3.25, 3, 2.25))

xlim=c(0.5,4.5)
ylim=c(0,1)

plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=ylim)


vioplot(phast.pc, h=0.02, add=T, axes=F, at=1, border="black", pchMed=21, colMed="black", colMed2="white", cex=0.95, col="indianred")
vioplot(phast.cited.more.lnc, h=0.02, add=T, axes=F, at=2, border="black", pchMed=21, colMed="black", colMed2="white", cex=0.95, col="dodgerblue4")
vioplot(phast.cited.once.lnc, h=0.02, add=T, axes=F, at=3, border="black", pchMed=21, colMed="black", colMed2="white", cex=0.95, col="steelblue")
vioplot(phast.other.lnc, h=0.02, add=T, axes=F, at=4, border="black", pchMed=21, colMed="black", colMed2="white", cex=0.95, col="lightblue")

axis(side=2, mgp=c(3,0.65,0))

axis(side=1, at=1:4, labels=rep("",4))

mtext("sequence conservation", side=2, line=2.5, cex=0.8)

mtext("e", font=2, line=0.95, at=-0.6)

##########################################################################

## legend plot


par(mar=c(0, 0, 0, 10.75))


plot(1, type="n", xlab="", ylab="", axes=F)

legend("topleft", legend=c("protein-coding"), fill=c("indianred"), xpd=NA, inset=c(0.01, -0.3), cex=1.2, bty="n")
legend("topleft", legend=c("cited lncRNAs (n>1)"), fill=c("dodgerblue4"), xpd=NA, inset=c(0.25, -0.3), cex=1.2, bty="n")
legend("topleft", legend=c("cited lncRNAs (n=1)"), fill=c("steelblue"), xpd=NA, inset=c(0.55, -0.3), cex=1.2, bty="n")
legend("topleft", legend=c("other lncRNAs"), fill=c("lightblue"), xpd=NA, inset=c(0.85, -0.3), cex=1.2, bty="n")


##########################################################################

dev.off()

##########################################################################
