##########################################################################

if(!("pathFigures"%in%objects())){
    pathRData="../../data_for_publication/RData/"
    pathFigures="../../data_for_publication/main_figures/"

    library(ade4)

    load=TRUE
    prepare=TRUE
}

##########################################################################

if(load){
    ## cited lncRNAs
    load(paste(pathRData, "data.PubMed.analysis.RData", sep=""))

    ## gene info
    load(paste(pathRData, "data.gene.info.RData", sep=""))

    ## sampleinfo
    load(paste(pathRData, "data.sample.info.RData",sep=""))

    ## TPM
    load(paste(pathRData, "data.expression.levels.MERiC.RData",sep=""))
    tpm.meric=tpm

    col.Edmondson=c("gray40", "gold", "darkorange", "indianred", "red")
    names(col.Edmondson)=c("liver", "ES1", "ES2", "ES3", "ES4")

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

    pca.list=list()

    samples=c(liver.samples$biopsyID, tumor.samples$tumor_biopsyID)
    eg=c(rep("liver", nrow(liver.samples)), paste("ES",tumor.samples$edmondson,sep=""))

    samples.liver=liver.samples$biopsyID
    samples.eg1=tumor.samples$tumor_biopsyID[which(tumor.samples$edmondson==1)]
    samples.eg2=tumor.samples$tumor_biopsyID[which(tumor.samples$edmondson==2)]
    samples.eg3=tumor.samples$tumor_biopsyID[which(tumor.samples$edmondson==3)]
    samples.eg4=tumor.samples$tumor_biopsyID[which(tumor.samples$edmondson==4)]

    for(genetype in c("lnc.cited.once", "lnc.cited.more", "other.lnc", "pc.cited.once", "pc.cited.more", "other.pc")){

        genes=get(genetype)
        genes=intersect(genes, rownames(tpm.meric))
        this.tpm=log2(tpm.meric[genes,]+1)

        this.pca=dudi.pca(t(this.tpm), center=T, scale=F, scannf=FALSE, nf=5)

        pca.list[[genetype]]=this.pca
    }

    prepare=FALSE
}

##########################################################################

pdf(paste(pathFigures, "SupplementaryFigure2.pdf", sep=""), width=6.85, height=6.75)

##########################################################################

## layout

m=matrix(rep(NA, 19*14), nrow=19)

n1=4
n2=3

for(i in 1:6){
    m[i,]=c(rep(1,n1), rep(2,n2), rep(7,n1), rep(8,n2))
}

for(i in 7:12){
    m[i,]=c(rep(3,n1), rep(4,n2), rep(9,n1), rep(10,n2))
}

for(i in 13:18){
    m[i,]=c(rep(5,n1), rep(6,n2), rep(11,n1), rep(12,n2))
}

m[19, ]= rep(13,14)

layout(m)

##########################################################################

genetypes=c("pc.cited.more", "pc.cited.once", "other.pc", "lnc.cited.more", "lnc.cited.once",  "other.lnc")

titles=c("protein-coding, cited >1", "protein-coding, cited 1", "protein-coding, not cited", "lncRNAs, cited >1", "lncRNAs, cited 1", "lncRNAs, not cited")
names(titles)=genetypes

labels=letters[1:6]
names(labels)=c("pc.cited.more", "lnc.cited.more", "pc.cited.once", "lnc.cited.once", "other.pc",   "other.lnc")

##########################################################################

for(genetype in genetypes){
    this.pca=pca.list[[genetype]]

    ## first factorial map

    min.coord.liver=min(as.numeric(this.pca$li[liver.samples$biopsyID,1]))
    mean.coord.tumor=mean(as.numeric(this.pca$li[tumor.samples$tumor_biopsyID,1]))

    sign=1

    if(min.coord.liver>mean.coord.tumor){
        sign=-1
    }

    par(mar=c(4,3.1,2.1,1.1))
    plot(sign*this.pca$li[samples,1], this.pca$li[samples,2], pch=20, col=col.Edmondson[eg], xlab="", ylab="", axes=F)

    axis(side=1, mgp=c(3,0.5,0))
    axis(side=2, mgp=c(3,0.65,0))

    box()

    explained=round(100*this.pca$eig/sum(this.pca$eig), digits=1)

    mtext(paste("PC1 (",explained[1],"% variance)", sep=""), side=1, line=2.15, cex=0.75)
    mtext(paste("PC2 (",explained[2],"% variance)", sep=""), side=2, line=2, cex=0.75)

    xlim=range(sign*this.pca$li[samples,1])
    mtext(titles[genetype], side=3, at=xlim[2]+diff(xlim)/10, cex=0.85, line=0.5)

    mtext(labels[genetype], side=3, line=1, cex=0.9, font=2, at=xlim[1]-diff(xlim)/3.5)


    ## coordinates on first axis
    par(mar=c(4,3.1,2.1,1.1))

    boxplot(sign*this.pca$li[samples.liver,1], sign*this.pca$li[samples.eg1,1], sign*this.pca$li[samples.eg2,1], sign*this.pca$li[samples.eg3,1], sign*this.pca$li[samples.eg4,1], col=col.Edmondson, pch=20, axes=F)

    axis(side=1, mgp=c(3,0.5,0), at=1:5, labels=rep("",5))
    axis(side=2, mgp=c(3,0.65,0))
    box()

    mtext("coordinate on PC1", side=2, line=2, cex=0.75)
    mtext("sample type", side=1, line=2.15, cex=0.75)

}

##########################################################################

## legend plot
par(mar=c(0.5,3.1,0.1,1.1))
plot(1, type="n", xlab="", ylab="", axes=F)

legend("topleft", legend="liver", fill=col.Edmondson[1], bty="n", cex=1.25)
legend("topleft", legend="tumor grade 1", fill=col.Edmondson[2], bty="n", cex=1.25, inset=c(0.1,0))
legend("topleft", legend="tumor grade 2", fill=col.Edmondson[3], bty="n", cex=1.25, inset=c(0.3,0))
legend("topleft", legend="tumor grade 3", fill=col.Edmondson[4], bty="n", cex=1.25, inset=c(0.5,0))
legend("topleft", legend="tumor grade 4", fill=col.Edmondson[5], bty="n", cex=1.25, inset=c(0.7,0))

##########################################################################

dev.off()

##########################################################################
