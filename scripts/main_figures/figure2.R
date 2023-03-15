##############################################################################

pathExpression="../../results/expression_estimation/"
pathDifferentialExpression="../../results/differential_expression/"
pathDiffExpTCGA="../../results/differential_expression_TCGA/Ensembl/"
pathAnnot="../../data/ensembl_annotations/"
pathSampleInfo="../../results/sample_info/"
pathPubMed="../../results/PubMed_search/"
pathResults="../../results/figures/main_figures/"

annot="AllTranscripts_Ensembl109_noMT_norRNA_nohaplo"
release=109

minFC=1.5
maxFDR=1e-3

library(scales)

########################################################################

if(prepare){

########################################################################
  
  geneinfo=read.table(paste(pathAnnot, "GeneInfo_Ensembl",release,".txt",sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
  rownames(geneinfo)=geneinfo$Gene.stable.ID
  
  pc=geneinfo$Gene.stable.ID[which(geneinfo$Gene.type=="protein_coding")]
  lnc=geneinfo$Gene.stable.ID[which(geneinfo$Gene.type=="lncRNA")]
  
########################################################################
  
  de.tl=read.table(paste(pathDifferentialExpression,annot, "/DifferentialExpression_Tumor_vs_Liver.txt", sep=""), h=T, stringsAsFactors=F, sep="\t")
  
  de.grades=read.table(paste(pathDifferentialExpression,annot, "/DifferentialExpression_EdmondsonGrade_34_vs_12.txt", sep=""), h=T, stringsAsFactors=F, sep="\t")

  de.tcga=read.table(paste(pathDiffExpTCGA, "DifferentialExpressionResults_Tumor_vs_NonTumor.txt",sep=""), h=T, stringsAsFactors=F, sep="\t")
  rownames(de.tcga)=de.tcga$GeneID
  de.tcga=de.tcga[,-1]
  
########################################################################

  up.genes.tl=rownames(de.tl)[which(de.tl$padj<maxFDR & de.tl$log2FoldChange>=log2(minFC))]
  down.genes.tl=rownames(de.tl)[which(de.tl$padj<maxFDR & de.tl$log2FoldChange<=log2(1/minFC))]
  tested.genes.tl=rownames(de.tl)[which(!is.na(de.tl$padj))]
  signif.tl=c(up.genes.tl, down.genes.tl)

  up.genes.grades=rownames(de.grades)[which(de.grades$padj<maxFDR & de.grades$log2FoldChange>=log2(minFC))]
  down.genes.grades=rownames(de.grades)[which(de.grades$padj<maxFDR & de.grades$log2FoldChange<=log2(1/minFC))]
  tested.genes.grades=rownames(de.grades)[which(!is.na(de.grades$padj))]
  signif.grades=c(up.genes.grades, down.genes.grades)
  
  up.genes.tcga=rownames(de.tcga)[which(de.tcga$padj<maxFDR & de.tcga$log2FoldChange>=log2(minFC))]
  down.genes.tcga=rownames(de.tcga)[which(de.tcga$padj<maxFDR & de.tcga$log2FoldChange<=log2(1/minFC))]
  tested.genes.tcga=rownames(de.tcga)[which(!is.na(de.tcga$padj))]
  signif.tcga=c(up.genes.tcga, down.genes.tcga)
  
########################################################################

  nb.cit=read.table(paste(pathPubMed, "number_of_citations_hepatocellular_carcinoma_Title.txt",sep=""),h=T, stringsAsFactors=F)

  cited.pc=nb.cit$GeneID[which(nb.cit$GeneType=="protein_coding")]
  cited.lnc=nb.cit$GeneID[which(nb.cit$GeneType=="lncRNA")]
  
  ########################################################################
  
  prepare=T
}

########################################################################
########################################################################

pdf(file=paste(pathResults, "Figure2.pdf",sep=""), width=4.6, height=4.25)

m=matrix(c(rep(NA,15*8)), nrow=15)

for(i in 1:9){
  m[i,]=c(rep(1, 4), rep(2, 4))
}


for(i in 10:15){
  m[i,]=c(rep(3, 4), rep(4, 4))
}
layout(m)

########################################################################

## tumor vs liver

par(mar=c(5.1, 4.1, 2.5, 1))
plot(de.tl[signif.tl, "log2FoldChange"], -log10(de.tl[signif.tl, "padj"]), pch=20, xlab="", ylab="", type="n", axes=F)
points(de.tl[intersect(signif.tl, pc),  "log2FoldChange"], -log10(de.tl[intersect(signif.tl, pc), "padj"]), pch=20, col=alpha("gray40", 0.5))
points(de.tl[intersect(signif.tl, lnc),  "log2FoldChange"], -log10(de.tl[intersect(signif.tl, lnc), "padj"]), pch=20, col=alpha("steelblue", 0.5))

axis(side=1, cex=1.1, mgp=c(3, 0.5, 0))
axis(side=2, cex=1.1, mgp=c(3, 0.75, 0))

mtext("log2(fold change)", side=1, line=2, cex=0.75)
mtext("-log10(FDR)", side=2, line=2.25, cex=0.75)

box()

abline(h=3, lty=3)
abline(v=c(log2(1/1.5),log2(1.5)), lty=3)

mtext("tumor vs. normal liver", side=3, cex=0.75, line=0.5)

mtext("A", side=3, at=-9.75, font=2, cex=0.95, line=0.75)

legend("topleft", c("protein-coding", "lncRNAs"), inset=0.015, bg="white", pch=20, col=c("gray40", "steelblue"), bty="y", box.col="white")

########################################################################

## tumor grades

par(mar=c(5.1, 4.1, 2.5, 1))
plot(de.grades[signif.grades, "log2FoldChange"], -log10(de.grades[signif.grades, "padj"]), pch=20, xlab="", ylab="", type="n", axes=F)
points(de.grades[intersect(signif.grades, pc),  "log2FoldChange"], -log10(de.grades[intersect(signif.grades, pc), "padj"]), pch=20, col=alpha("gray40", 0.5))
points(de.grades[intersect(signif.grades, lnc),  "log2FoldChange"], -log10(de.grades[intersect(signif.grades, lnc), "padj"]), pch=20, col=alpha("steelblue", 0.5))

axis(side=1, cex=1.1, mgp=c(3, 0.5, 0))
axis(side=2, cex=1.1, mgp=c(3, 0.75, 0))

mtext("log2(fold change)", side=1, line=2, cex=0.75)
mtext("-log10(FDR)", side=2, line=2.25, cex=0.75)

box()

abline(h=3, lty=3)
abline(v=c(log2(1/1.5),log2(1.5)), lty=3)

mtext("B", side=3, at=-7.15, font=2, cex=0.95, line=0.75)

mtext("tumor grades", side=3, cex=0.75, line=0.5)

########################################################################

## % DE tl

par(mar=c(4.5, 4.1, 0.5, 1))

values.tl=100*c(length(intersect(pc, up.genes.tl))/length(intersect(pc, tested.genes.tl)), length(intersect(pc, down.genes.tl))/length(intersect(pc, tested.genes.tl)), length(intersect(lnc, up.genes.tl))/length(intersect(lnc, tested.genes.tl)), length(intersect(lnc, down.genes.tl))/length(intersect(lnc, tested.genes.tl)), length(intersect(cited.pc, up.genes.tl))/length(intersect(cited.pc, tested.genes.tl)), length(intersect(cited.pc, down.genes.tl))/length(intersect(cited.pc, tested.genes.tl)), length(intersect(cited.lnc, up.genes.tl))/length(intersect(cited.lnc, tested.genes.tl)), length(intersect(cited.lnc, down.genes.tl))/length(intersect(cited.lnc, tested.genes.tl)))

b=barplot(values.tl, space=c(0.1, 0.2, 0.5, 0.2, 0.75, 0.2, 0.5, 0.2), col=rep(c("gray40","gray40", "steelblue", "steelblue"),2), axes=F, ylim=c(0,15))
axis(side=2, mgp=c(3, 0.5, 0), cex.axis=1, at=c(0,5,10, 15))
axis(side=1, at=b, labels=rep("", 8))

mtext("% DE genes", side=2, line=2.25, cex=0.75)

mtext(rep(c("up", "down"),4), at=b, side=1, las=2, cex=0.6, line=0.75)

abline(lty=3, v=mean(b[4:5]))

##legend("topleft", c("protein-coding", "lncRNAs"), inset=c(0.015, -0.25), border="black", fill=c("gray40", "steelblue"), xpd=NA, bty="y", box.col="white", bg="white")

mtext("all genes", side=1, at=mean(b[2:3]), line=3.25, cex=0.7)
mtext("HCC-associated", side=1, at=mean(b[6:7]), line=3.25, cex=0.7)

mtext("C", side=3, at=-2.9, font=2, cex=0.95, line=0.85)


########################################################################

## % DE tl

par(mar=c(4.5, 4.1, 0.5, 1))

values.grades=100*c(length(intersect(pc, up.genes.grades))/length(intersect(pc, tested.genes.grades)), length(intersect(pc, down.genes.grades))/length(intersect(pc, tested.genes.grades)), length(intersect(lnc, up.genes.grades))/length(intersect(lnc, tested.genes.grades)), length(intersect(lnc, down.genes.grades))/length(intersect(lnc, tested.genes.grades)), length(intersect(cited.pc, up.genes.grades))/length(intersect(cited.pc, tested.genes.grades)), length(intersect(cited.pc, down.genes.grades))/length(intersect(cited.pc, tested.genes.grades)), length(intersect(cited.lnc, up.genes.grades))/length(intersect(cited.lnc, tested.genes.grades)), length(intersect(cited.lnc, down.genes.grades))/length(intersect(cited.lnc, tested.genes.grades)))

b=barplot(values.grades, space=c(0.1, 0.2, 0.5, 0.2, 0.75, 0.2, 0.5, 0.2), col=rep(c("gray40","gray40", "steelblue", "steelblue"),2), axes=F, ylim=c(0,15))
axis(side=2, mgp=c(3, 0.5, 0), cex.axis=1, at=c(0,5,10, 15))
axis(side=1, at=b, labels=rep("", 8))

mtext("% DE genes", side=2, line=2.25, cex=0.75)

mtext(rep(c("up", "down"),4), at=b, side=1, las=2, cex=0.6, line=0.75)

abline(lty=3, v=mean(b[4:5]))


mtext("all genes", side=1, at=mean(b[2:3]), line=3.25, cex=0.7)
mtext("HCC-associated", side=1, at=mean(b[6:7]), line=3.25, cex=0.7)

mtext("D", side=3, at=-2.9, font=2, cex=0.95, line=0.85)

########################################################################

dev.off()

########################################################################
