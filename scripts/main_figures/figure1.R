##############################################################################

pathExpression="../../results/expression_estimation/"
pathDifferentialExpression="../../results/differential_expression/"
pathAnnot="../../data/ensembl_annotations/"
pathSampleInfo="../../results/sample_info/"
pathResults="../../results/figures/main_figures/"

annot="AllTranscripts_Ensembl109_noMT_norRNA_nohaplo"

set.seed(19)

library(ade4)

########################################################################

if(prepare){

  tpm=read.table(paste(pathExpression, annot, "/AllSamples_KallistoRawTPM.txt",sep=""), h=T, stringsAsFactors=F)
  
  ########################################################################
  
  sampleinfo=read.table(paste(pathSampleInfo, "AnalyzedSamples.txt", sep=""), h=T, stringsAsFactors=F,sep="\t")
  
  col.samples=rep("gray40", nrow(sampleinfo))
  col.samples[which(sampleinfo$TissueType=="Tumor" & sampleinfo$EdmondsonGrade==1)]="yellow"
  col.samples[which(sampleinfo$TissueType=="Tumor" & sampleinfo$EdmondsonGrade==2)]="orange"
  col.samples[which(sampleinfo$TissueType=="Tumor" & sampleinfo$EdmondsonGrade==3)]="darkorange"
  col.samples[which(sampleinfo$TissueType=="Tumor" & sampleinfo$EdmondsonGrade==4)]="red"
  
  names(col.samples)=sampleinfo$BiopsyID
  
########################################################################
  
  geneinfo=read.table(paste(pathAnnot, "GeneInfo_Ensembl109.txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
  
  pc=geneinfo$Gene.stable.ID[which(geneinfo$Gene.type=="protein_coding")]
  
  lnc=geneinfo$Gene.stable.ID[which(geneinfo$Gene.type=="lncRNA")]
  
  pseudo=geneinfo$Gene.stable.ID[which(geneinfo$Gene.type%in%c("transcribed_unitary_pseudogene", "transcribed_unprocessed_pseudogene", "unitary_pseudogene", "unprocessed_pseudogene"))]
  
 ########################################################################
  
  tpm.pc=tpm[intersect(pc, rownames(tpm)),]
  
  tpm.lnc=tpm[intersect(lnc, rownames(tpm)),]
  
  tpm.pseudo=tpm[intersect(pseudo, rownames(tpm)),]
  
 ########################################################################
  
  pca.pc=dudi.pca(t(log2(tpm.pc+1)), scale=F, center=T, scannf=F, nf=5)
  
  pca.lnc=dudi.pca(t(log2(tpm.lnc+1)), scale=F, center=T, scannf=F, nf=5)
  
  pca.pseudo=dudi.pca(t(log2(tpm.pseudo+1)), scale=F, center=T, scannf=F, nf=5)

  explained.pc=round(100*pca.pc$eig/sum(pca.pc$eig))

  explained.lnc=round(100*pca.lnc$eig/sum(pca.lnc$eig))
  
  explained.pseudo=round(100*pca.pseudo$eig/sum(pca.pseudo$eig))
  

  ########################################################################

  prepare=F
}

########################################################################
#######################################################################

pdf(file=paste(pathResults, "Figure1.pdf",sep=""),width=10,height=4)

m=matrix(c(rep(NA,10*9)), nrow=10)

for(i in 1:10){
  m[i,]=c(rep(1, 3), rep(2, 3), rep(3, 3))
}
layout(m)

par(mar=c(4.1,4.1,2.1,1.1))

plot(-pca.pc$li[,1], pca.pc$li[,2], pch=20, col=col.samples[rownames(pca.pc$li)], xlab="", ylab="", axes=F, cex=1.25)
axis(side=1, mgp=c(3, 0.75, 0), cex.axis=1.2)
axis(side=2, mgp=c(3, 0.5, 0), cex.axis=1.2)
box()
mtext("protein-coding genes", side=3, line=0.25)
mtext(paste("PC1 (", explained.pc[1],"% variance)",sep=""), side=1, line=2.25)
mtext(paste("PC2 (", explained.pc[2],"% variance)",sep=""), side=2, line=2.25)

plot(pca.lnc$li[,1], pca.lnc$li[,2], pch=20, col=col.samples[rownames(pca.lnc$li)], ylab="", xlab="",axes=F, cex=1.25)
axis(side=1, mgp=c(3, 0.75, 0), cex.axis=1.2)
axis(side=2, mgp=c(3, 0.5, 0), cex.axis=1.2)
box()
mtext("lncRNAs", side=3, line=0.25)
mtext(paste("PC1 (", explained.lnc[1],"% variance)",sep=""), side=1, line=2.25)
mtext(paste("PC2 (", explained.lnc[2],"% variance)",sep=""), side=2, line=2.25)


plot(pca.pseudo$li[,1], -pca.pseudo$li[,2], pch=20, col=col.samples[rownames(pca.pseudo$li)], ylab="", xlab="", axes=F, cex=1.25)
axis(side=1, mgp=c(3, 0.75, 0), cex.axis=1.2)
axis(side=2, mgp=c(3, 0.5, 0), cex.axis=1.2)
box()
mtext("pseudogenes", side=3, line=0.25)
mtext(paste("PC1 (", explained.pseudo[1],"% variance)",sep=""), side=1, line=2.25)
mtext(paste("PC2 (", explained.pseudo[2],"% variance)",sep=""), side=2, line=2.25)



dev.off()

########################################################################
