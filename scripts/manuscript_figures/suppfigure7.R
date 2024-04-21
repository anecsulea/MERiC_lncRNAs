#############################################################################

if(!("pathFigures"%in%objects())){
    pathRData="../../data_for_publication/RData/"
    pathFigures="../../data_for_publication/main_figures/"

    library(vioplot)

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

     ## sample info
    load(paste(pathRData, "data.sample.info.TCGA.RData",sep=""))

    ## TPM
    load(paste(pathRData, "data.expression.levels.TCGA.RData",sep=""))
    tpm.tcga=tpm


    maxFDR=0.05
    minFC=1.5

    load=FALSE
}

#############################################################################

if(prepare){
  ## for each gene, how many articles mention it

    highly.cited.lnc=names(nb.citations.lnc)[which(nb.citations.lnc>=10)]

    highly.cited.lnc=highly.cited.lnc[order(nb.citations.lnc[highly.cited.lnc], decreasing=T)]

    highly.cited.lnc=setdiff(highly.cited.lnc, "ENSG00000235300")
    ## synonym of this lncRNA is BTR, which refers to BCAA/tyrosine ratio in the HCC  literature

    ## DE table

    de.table=matrix(rep(NA, length(highly.cited.lnc)*1), nrow=length(highly.cited.lnc))

    rownames(de.table)=highly.cited.lnc
    colnames(de.table)=c("tnt.tcga")

    for(test in colnames(de.table)){
        this.de=get(paste("diffexp", test, sep="."))
        tested=rownames(this.de)
        signif.up=rownames(this.de)[which(this.de$padj<maxFDR & this.de$log2FoldChange>log2(minFC))]
        signif.down=rownames(this.de)[which(this.de$padj<maxFDR & this.de$log2FoldChange<log2(1/minFC) )]

        relax.up=rownames(this.de)[which(this.de$padj<maxFDR & this.de$log2FoldChange>0)]
        relax.down=rownames(this.de)[which(this.de$padj<maxFDR & this.de$log2FoldChange<0 )]

        relax.up=setdiff(relax.up, signif.up)
        relax.down=setdiff(relax.down, signif.down)

        if(any(highly.cited.lnc%in%tested)){
            de.table[intersect(highly.cited.lnc, tested), test]=0
        }

        if(any(highly.cited.lnc%in%signif.up)){
            de.table[intersect(highly.cited.lnc, signif.up), test]=1
        }
        if(any(highly.cited.lnc%in%relax.up)){
            de.table[intersect(highly.cited.lnc, relax.up), test]=0.5
        }

        if(any(highly.cited.lnc%in%signif.down)){
            de.table[intersect(highly.cited.lnc, signif.down), test]=-1
        }

        if(any(highly.cited.lnc%in%relax.down)){
            de.table[intersect(highly.cited.lnc, relax.down), test]=-0.5
        }
    }

    ## compute expression difference between paired biopsies

    tumor.samples=rownames(sampleinfo)[which(sampleinfo$sample_type=="Tumor")]
    names(tumor.samples)=sampleinfo$case_id[which(sampleinfo$sample_type=="Tumor")]

    nontumor.samples=rownames(sampleinfo)[which(sampleinfo$sample_type=="Non-Tumor")]
    names(nontumor.samples)=sampleinfo$case_id[which(sampleinfo$sample_type=="Non-Tumor")]

    print(all(names(tumor.samples)==names(nontumor.samples)))

    tpm.tumor=tpm.tcga[, tumor.samples]
    tpm.nontumor=tpm.tcga[, nontumor.samples]

    unique.patients=names(tumor.samples)

    prepare=FALSE
}

#############################################################################

## 1 column width 87 mm = 3.34 in
## 2 columns width 180 mm = 7.04 in
## max height: 9 in, including legend

#############################################################################

pdf(file=paste(pathFigures, "SupplementaryFigure7.pdf", sep=""), width=3.85, height=7.5)

m=matrix(rep(NA, 28*6), nrow=28)

for(i in seq(from=1, by=1, length.out=length(highly.cited.lnc))){
    j=2*i-1
    m[i,]=c(rep(j,5), rep(j+1,1))
}

## legend

for(i in 27:28){
  M=max(m, na.rm=T)
  m[i,]=rep(M+1, 6)
}

layout(m)

#############################################################################

## selected lncRNAs

for(id in highly.cited.lnc){
    ## gene name
    name=geneinfo[id,"Name"]

    print(name)

    ## gene expression values for sample categories
    this.exp.tumor=as.numeric(tpm.tumor[id, ])
    this.exp.nontumor=as.numeric(tpm.nontumor[id, ])

    this.meanexp=(this.exp.tumor+this.exp.nontumor)/2

    this.diffexp=(this.exp.tumor-this.exp.nontumor)/this.meanexp

    ## difference between tumor and liver

    par(mar=c(0.5, 6.5, 0.1, 0.15))

    ## xlim=range(c(this.exp.nontumor, this.exp.tumor))

    xlim=c(-2,2)
    ylim=c(0.25, 1.75)
    plot(1, type="n", xlim=xlim, ylim=ylim, axes=F, xlab="", ylab="")

    vioplot(this.diffexp, h=0.25, add=T, axes=F, at=1, border="black", pchMed=21, colMed="black", colMed2="white", cex=0.95, col="gray60", horizontal=T)

    ## axis if last plot
    if(id==highly.cited.lnc[length(highly.cited.lnc)]){
        axis(side=1, at=seq(from=-2, to=2, by=2), mgp=c(3, 0.5, 0), cex=0.7)
        axis(side=1, at=seq(from=-1, to=1, by=2), mgp=c(3, 0.5, 0), cex=0.7)

        segments(0, 0, 0, 62, lty=2, col="darkred", xpd=NA)

        mtext("TPM normalized difference, tumor vs. adjacent tissue", side=1, at=0, line=2.5, cex=0.75)
    }

    mtext(name, side=2, las=2, cex=0.6, line=6.2, font=3, adj=0)

    abline(h=0.2, lty=3, col="gray40", xpd=NA)

    ## differential expression

    par(mar=c(0.5,1.65,0.1,0))
    xlim=c(0,1)
    ylim=c(0,1)
    plot(1, type="n", xlab="", ylab="",axes=F, xlim=xlim, ylim=ylim, xaxs="i", yaxs="i")

    ## tumor vs normal liver

    x2=0.5

    ## tumor vs non-tumor, paired samples

    if(de.table[id,"tnt.tcga"]==1){
        arrows(x2, 0.1, x2, 0.9, length=0.035, lwd=1.15, col="black")
    }
    if(de.table[id,"tnt.tcga"]==0.5){
        arrows(x2, 0.1, x2, 0.9, length=0.035, lwd=1.15, col="gray50")
    }

    if(de.table[id,"tnt.tcga"]==-1){
        arrows(x2, 0.9, x2, 0.1, length=0.035, lwd=1.15, col="black")
    }
    if(de.table[id,"tnt.tcga"]==-0.5){
        arrows(x2, 0.9, x2, 0.1, length=0.035, lwd=1.15, col="gray50")
    }

}

#############################################################################

dev.off()

#############################################################################
