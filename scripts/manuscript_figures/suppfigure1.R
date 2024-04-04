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
    other.lnc=setdiff(lnc, lnc.cited.all)

    pc.cited.once=names(nb.citations.pc)[which(nb.citations.pc==1)]
    pc.cited.more=names(nb.citations.pc)[which(nb.citations.pc>1)]
    pc.cited.all=names(nb.citations.pc)
    other.pc=setdiff(pc, pc.cited.all)

    ## average expression level across various samples

    meantpm.tumor.tcga=apply(tpm.tcga[, tumor.samples.tcga$id], 1, mean)
    meantpm.nontumor.tcga=apply(tpm.tcga[, nontumor.samples.tcga$id], 1, mean)

    prepare=FALSE
}

##########################################################################

pdf(paste(pathFigures, "SupplementaryFigure1.pdf", sep=""), width=3.85, height=3)

##########################################################################

## layout

m=matrix(rep(NA, 1*10), nrow=10)

for(i in 1:10){
    m[i,]=c(rep(1, 10), rep(2,10))
}


layout(m)

##########################################################################

genetypes=c("pc.cited.more", "pc.cited.once", "other.pc", "lnc.cited.more", "lnc.cited.once",  "other.lnc")

xpos.genetypes=c(1, 3.5, 6, 2,  4.5, 7)
names(xpos.genetypes)=genetypes

col.genetypes=rep(c("indianred", "steelblue"), each=3)
names(col.genetypes)=genetypes

##########################################################################

## expression levels in TCGA, tumor samples

par(mar=c(2.5, 3.75, 2.5, 1.75))

xlim=c(0.5,8)
ylim=c(0,20)

plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=ylim)

## expression in tumor samples

for(type in genetypes){
    this.genes=get(type)
    this.xpos=xpos.genetypes[type]
    this.col=col.genetypes[type]

    vioplot(log2(meantpm.tumor.tcga[this.genes]+1), h=1, add=T, axes=F, at=this.xpos, border="black", pchMed=21, colMed="black", colMed2="white", cex=0.95, col=this.col)
}

abline(v=mean(xpos.genetypes[c("lnc.cited.more", "pc.cited.once")]), lty=3)
abline(v=mean(xpos.genetypes[c("lnc.cited.once", "other.pc")]), lty=3)

mtext("tumor", side=3, line=-1, at=mean(xpos.genetypes), cex=0.75)

axis(side=2, mgp=c(3,0.65,0))
mtext("mean expression level (log2 TPM)", side=2, line=2.5, cex=0.75)

mtext("cited", side=1, at=mean(xpos.genetypes[c("pc.cited.more", "lnc.cited.more")]),line=0.5, cex=0.75)
mtext(">1 articles", side=1, at=mean(xpos.genetypes[c("pc.cited.more", "lnc.cited.more")]),line=1.5, cex=0.75)

mtext("cited", side=1, at=mean(xpos.genetypes[c("pc.cited.once", "lnc.cited.once")]),line=0.5, cex=0.75)
mtext("1 article", side=1, at=mean(xpos.genetypes[c("pc.cited.once", "lnc.cited.once")]),line=1.5, cex=0.75)

mtext("not cited", side=1, at=mean(xpos.genetypes[c("other.pc", "other.lnc")]),line=1, cex=0.75)

axis(side=1, at=xpos.genetypes, labels=rep("", length(xpos.genetypes)), mgp=c(3,0.5,0))

legend("topright", legend=c("protein-coding", "lncRNAs"), fill=c("indianred", "steelblue"), xpd=NA, inset=c(-0.05, -0.05), cex=1.1, bty="n")

mtext("a", font=2, line=0.95, at=-1)

##########################################################################

## expression levels in TCGA, adjacent tissue samples

par(mar=c(2.5, 4.75, 2.5, 0.75))

xlim=c(0.5,8)
ylim=c(0,20)

plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=ylim)

## expression in tumor samples

for(type in genetypes){
    this.genes=get(type)
    this.xpos=xpos.genetypes[type]
    this.col=col.genetypes[type]

    vioplot(log2(meantpm.nontumor.tcga[this.genes]+1), h=1, add=T, axes=F, at=this.xpos, border="black", pchMed=21, colMed="black", colMed2="white", cex=0.95, col=this.col)
}

abline(v=mean(xpos.genetypes[c("lnc.cited.more", "pc.cited.once")]), lty=3)
abline(v=mean(xpos.genetypes[c("lnc.cited.once", "other.pc")]), lty=3)

mtext("adjacent tissue", side=3, line=-1, at=mean(xpos.genetypes), cex=0.75)

axis(side=2, mgp=c(3,0.65,0))
mtext("mean expression level (log2 TPM)", side=2, line=2.5, cex=0.75)

mtext("cited", side=1, at=mean(xpos.genetypes[c("pc.cited.more", "lnc.cited.more")]),line=0.5, cex=0.75)
mtext(">1 articles", side=1, at=mean(xpos.genetypes[c("pc.cited.more", "lnc.cited.more")]),line=1.5, cex=0.75)

mtext("cited", side=1, at=mean(xpos.genetypes[c("pc.cited.once", "lnc.cited.once")]),line=0.5, cex=0.75)
mtext("1 article", side=1, at=mean(xpos.genetypes[c("pc.cited.once", "lnc.cited.once")]),line=1.5, cex=0.75)

mtext("not cited", side=1, at=mean(xpos.genetypes[c("other.pc", "other.lnc")]),line=1, cex=0.75)

axis(side=1, at=xpos.genetypes, labels=rep("", length(xpos.genetypes)), mgp=c(3,0.5,0))

mtext("b", font=2, line=0.95, at=-1)

##########################################################################

dev.off()

##########################################################################
