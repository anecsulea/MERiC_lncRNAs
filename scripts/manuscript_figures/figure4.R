#############################################################################

if(!("pathFigures"%in%objects())){
    pathRData="../../data_for_publication/RData/"
    pathFigures="../../data_for_publication/figures/"

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
    load(paste(pathRData, "data.sample.info.MERiC.RData",sep=""))

    ## TPM
    load(paste(pathRData, "data.expression.levels.MERiC.RData",sep=""))
    tpm.meric=tpm

    col.Edmondson=c("gold", "darkorange", "indianred", "red")

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

    ## sample info, only for tumor samples

    rownames(tumor.samples)=as.character(tumor.samples$tumor_biopsyID)

    tumor.samples$NAFLD=rep("no", nrow(tumor.samples))
    tumor.samples$NAFLD[grep("NAFLD", tumor.samples$underlying_liver_disease)]="yes"

    tumor.samples$ALD=rep("no", nrow(tumor.samples))
    tumor.samples$ALD[grep("ALD", tumor.samples$underlying_liver_disease)]="yes"

    tumor.samples$HepB=rep("no", nrow(tumor.samples))
    tumor.samples$HepB[grep("Hepatitis B", tumor.samples$underlying_liver_disease)]="yes"

    tumor.samples$HepC=rep("no", nrow(tumor.samples))
    tumor.samples$HepC[grep("Hepatitis C", tumor.samples$underlying_liver_disease)]="yes"

    ## we order them by all the factors
    tumor.samples=tumor.samples[order(tumor.samples$sex),]
    tumor.samples=tumor.samples[order(tumor.samples$NAFLD),]
    tumor.samples=tumor.samples[order(tumor.samples$ALD),]
    tumor.samples=tumor.samples[order(tumor.samples$HepB),]
    tumor.samples=tumor.samples[order(tumor.samples$HepC),]
    tumor.samples=tumor.samples[order(tumor.samples$cirrhosis),]
    tumor.samples=tumor.samples[order(tumor.samples$edmondson),]

    ## DE table

    de.table=matrix(rep(NA, length(highly.cited.lnc)*2), nrow=length(highly.cited.lnc))

    rownames(de.table)=highly.cited.lnc
    colnames(de.table)=c("tnt.meric", "grades")

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

    tpm.tumor=tpm.meric[, tumor.samples$tumor_biopsyID]
    meantpm.tumor.patient=t(apply(tpm.tumor, 1, function(x) tapply(as.numeric(x), as.factor(tumor.samples$Patient_ID), mean)))


    tpm.nontumor=tpm.meric[, nontumor.samples$biopsyID]
    meantpm.nontumor.patient=t(apply(tpm.nontumor,1, function(x) tapply(as.numeric(x), as.factor(nontumor.samples$Patient_ID), mean)))

    unique.patients=unique(tumor.samples$Patient_ID)

    prepare=FALSE
}

#############################################################################

## 1 column width 87 mm = 3.34 in
## 2 columns width 180 mm = 7.04 in
## max height: 9 in, including legend

#############################################################################

pdf(file=paste(pathFigures, "Figure4.pdf", sep=""), width=6.85, height=9)

m=matrix(rep(NA, 62*19), nrow=62)

## first line: spacer
for(i in 1){
  m[i,]=rep(1, 19)
}

## individual lines:  Edmondson grade, cirrhosis, hepC, hepB, ALD, NAFLD, sex

for(i in 2:8){
  m[i,]=c(rep(2*(i-1), 6),  rep(2*(i-1)+1, 13))
}

## spacer plot

for(i in 9:10){
  m[i,]=rep(16, 19)
}

## expression of highly cited lncRNAs


## a*11+b=17
## a*13+b=20

a=(20-17)/2
b=17-11*a

for(i in seq(from=11, by=2, length.out=length(highly.cited.lnc))){
    j=a*i+b

  m[i,]=c(rep(j,5), rep(j+1,1),  rep(j+2,13))
  m[i+1,]=c(rep(j,5), rep(j+1,1),  rep(j+2,13))
}


## legend

for(i in 61:62){
  M=max(m, na.rm=T)
  m[i,]=rep(M+1, 19)
}

layout(m)

#############################################################################

## legend

par(mar=c(0,0,0,0))
plot(1, type="n", xlab="", ylab="", axes=F, xlim=c(0,1))

mtext("a", side=3, at=-0.025, line=-1.2, cex=0.95, font=2)

#############################################################################

labels=c("Edmondson grade",  "Cirrhosis", "Hepatitis C", "Hepatitis B", "ALD", "NAFLD", "Sex")
names(labels)=c("edmondson",  "cirrhosis", "HepC", "HepB", "ALD", "NAFLD","sex")

for(column in (names(labels))){
    if(column=="edmondson"){
        col.factor=col.Edmondson
        this.colors=col.Edmondson[tumor.samples[,"edmondson"]]
        levels=1:4
    } else{

        col.factor=gray.colors(n=length(levels(as.factor(tumor.samples[,column]))), rev=TRUE)
        names(col.factor)=levels(as.factor(tumor.samples[,column]))

        this.colors=col.factor[as.character(tumor.samples[,column])]

        levels=names(col.factor)
    }

  label=labels[column]

  ## legend for this plot

  par(mar=c(0,0,0,0))
  plot(1, type="n", xlab="", ylab="", axes=F, ylim=c(0,1), xlim=c(0, 1), xaxs="i", yaxs="i")

  ## mtext(label, side=3, cex=0.5, line=0, at=0)

  if(column=="BiopsyAge"){
    text(label, x=0.1, y=0.5, adj=c(0, 0.5), cex=0.95, xpd=NA)
    legend("topleft", legend=levels[1], fill=col.age[1], border=col.age[1], bty="n", horiz=T, xpd=NA, cex=0.85, inset=c(0.2, -0.3))
    legend("topleft", legend=levels[2:4], fill=col.age[2:4], border=col.age[2:4], bty="n", horiz=T, xpd=NA, cex=0.85, inset=c(0.4, -0.3))
    legend("topleft", legend=levels[5], fill=col.age[5], border=col.age[5], bty="n", horiz=T, xpd=NA, cex=0.85, inset=c(0.2, 0.55))
  } else{
    if(column=="edmondson"){
      legend("topleft", legend=levels, fill=col.factor, border=col.factor, bty="n", horiz=T, xpd=NA, cex=0.85, inset=c(0.55, -0.32))
      text(label, x=0.1, y=0.55, adj=c(0, 0.5), cex=0.95, xpd=NA)
    } else{
      legend("topleft", legend=levels[2:1], fill=col.factor[2:1], border=col.factor[2:1], bty="n", horiz=T, xpd=NA, cex=0.85, inset=c(0.55, -0.32))
      text(label, x=0.1, y=0.55, adj=c(0, 0.5), cex=0.95, xpd=NA)
     }
  }

  ## actual plot

  par(mar=c(0.1, 2.8, 0.15, 3))
  image(cbind(1:dim(tumor.samples)[1]), col=this.colors, axes=F, xaxs="i", yaxs="i")

}

#############################################################################

## spacer plot

plot(1, type="n", xlab="", ylab="", axes=F, xlim=c(0,1))
mtext("b", side=3, at=-0.087, line=-1.25, cex=0.95, font=2)
mtext("c", side=3, at=0.25, line=-1.25, cex=0.95, font=2)
mtext("d", side=3, at=0.31, line=-1.25, cex=0.95, font=2)
mtext("e", side=3, at=1.055, line=-1.25, cex=0.95, font=2)

#############################################################################

## selected lncRNAs

for(id in highly.cited.lnc){
    ## gene name
    name=geneinfo[id,"Name"]

    print(name)

    ## gene expression values for sample categories
    this.exp.tumor=as.numeric(meantpm.tumor.patient[id, unique.patients])
    this.exp.nontumor=as.numeric(meantpm.nontumor.patient[id, unique.patients])

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

        segments(0, 0, 0, 55, lty=2, col="darkred", xpd=NA)
    }

    mtext(name, side=2, las=2, cex=0.6, line=6.2, font=3, adj=0)

    abline(h=0.2, lty=3, col="gray40", xpd=NA)
   ## segments(-10, 0.25, xlim[2]+diff(xlim)/10, 0.25, lty=3, xpd=NA)

    ## differential expression

    par(mar=c(0.5,1.65,0.1,0))
    xlim=c(0,1)
    ylim=c(0,1)
    plot(1, type="n", xlab="", ylab="",axes=F, xlim=xlim, ylim=ylim, xaxs="i", yaxs="i")

    ## tumor vs normal liver

    x2=0.5

    ## tumor vs non-tumor, paired samples

    if(de.table[id,"tnt.meric"]==1){
        arrows(x2, 0.1, x2, 0.9, length=0.035, lwd=1.15, col="black")
    }
    if(de.table[id,"tnt.meric"]==0.5){
        arrows(x2, 0.1, x2, 0.9, length=0.035, lwd=1.15, col="gray50")
    }

    if(de.table[id,"tnt.meric"]==-1){
        arrows(x2, 0.9, x2, 0.1, length=0.035, lwd=1.15, col="black")
    }
    if(de.table[id,"tnt.meric"]==-0.5){
        arrows(x2, 0.9, x2, 0.1, length=0.035, lwd=1.15, col="gray50")
    }

    ## plot for tumors only

    this.exp.tumor=as.numeric(tpm.meric[id, tumor.samples$tumor_biopsyID])

    ylim=c(0, max(c(log2(this.exp.tumor+1),1)))

    xlim=c(0.25, length(this.exp.tumor)+0.5)
    xwidth=0.25

    par(mar=c(0.5, 2.8, 0.1, 3))
    plot(1, type="n", xlab="", ylab="",axes=F, xlim=xlim, ylim=ylim, xaxs="i", yaxs="i")

    xpos=1:length(this.exp.tumor)
    this.values=log2(this.exp.tumor+1)
    this.colors=col.Edmondson[tumor.samples$edmondson]
    rect(xpos-xwidth, 0, xpos+xwidth, this.values, col=this.colors, border=this.colors)

    ## abline(v=max(xpos), xpd=NA, lwd=0.2) ## to check alignment
    ## abline(v=min(xpos), xpd=NA) ## to check alignment

    maxval=ylim[2]
    par(tck=-0.25)
    axis(side=1, at=range(xlim)+c(-1,1), labels=rep("",2))
    axis(side=2, at=c(0, ylim[2]), labels=rep("",2), cex.axis=0.5, mgp=c(3, 0.5, 0))
    mtext(c(0, round(ylim[2], digits=1)), at=c(0, ylim[2]), las=2, cex=0.5, line=0.65, adj=1, side=2)

    ## text(name, x=xlim[2], adj=1, y=ylim[2]+diff(ylim)/5, cex=0.75, xpd=NA)

    ## DE, Edmondson grade

    x3=xlim[2]+diff(xlim)/20
    y1=diff(ylim)/10
    y2=ylim[2]-diff(ylim)/10

    if(de.table[id,"grades"]==1){
        arrows(x3, y1, x3, y2, length=0.035, lwd=1.15, xpd=NA, col="black")
    }
    if(de.table[id,"grades"]==0.5){
        arrows(x3, y1, x3, y2, length=0.035, lwd=1.15, xpd=NA, col="gray50")
    }

    if(de.table[id,"grades"]==-1){
        arrows(x3, y2, x3, y1, length=0.035, lwd=1.15, xpd=NA, col="black")
    }
    if(de.table[id,"grades"]==-0.5){
        arrows(x3, y2, x3, y1, length=0.035, lwd=1.15, xpd=NA, col="gray50")
    }
}

#############################################################################

par(mar=c(0, 0, 0, 0))
plot(1, type="n", xlab="", ylab="", axes=F, xlim=c(0,1), ylim=c(0,1))

text("1 & 2", x=1.01, y=1.1, cex=0.9, xpd=NA)
text("vs.", x=1.01, y=0.38, cex=0.9, xpd=NA, font=3)
text("3 & 4", x=1.01, y=-0.36, cex=0.9, xpd=NA)

text("TPM", x=0.05, y=0.45, cex=0.9, xpd=NA)
text("normalized difference", x=0.05, y=-0.35, cex=0.9, xpd=NA)


text("tumor", x=0.288, y=0.45, cex=0.9, xpd=NA)
text("vs.", x=0.288, y=-0.2, cex=0.9, xpd=NA, font=3)
text("adjacent tissue", x=0.288, y=-0.8, cex=0.9, xpd=NA)

legend("topleft", legend=as.character(1:4), fill=col.Edmondson, border=col.Edmondson, bty="n", cex=0.85, inset=c(0.77, 0), horiz=T, xpd=NA)



text("log2(TPM+1)", x=0.37, y=0.12, xpd=NA, cex=0.9)
text("Edmondson grade", x=0.72, y=0.3, cex=0.9, xpd=NA)

#############################################################################

dev.off()

#############################################################################
