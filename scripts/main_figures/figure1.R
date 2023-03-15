##############################################################################

pathExpression="../../results/expression_estimation/"
pathDifferentialExpression="../../results/differential_expression/"
pathAnnot="../../data/ensembl_annotations/"
pathSampleInfo="../../results/sample_info/"
pathResults="../../results/figures/main_figures/"

annot="AllTranscripts_Ensembl109_noMT_norRNA_nohaplo"

set.seed(19)

library(ade4)
library(vioplot)

########################################################################

if(prepare){

  ########################################################################
  
  sampleinfo=read.table(paste(pathSampleInfo, "AnalyzedSamples.txt", sep=""), h=T, stringsAsFactors=F,sep="\t")
  
  col.samples=rep("gray40", nrow(sampleinfo))
  col.samples[which(sampleinfo$TissueType=="Tumor" & sampleinfo$EdmondsonGrade==1)]="gold"
  col.samples[which(sampleinfo$TissueType=="Tumor" & sampleinfo$EdmondsonGrade==2)]="orange"
  col.samples[which(sampleinfo$TissueType=="Tumor" & sampleinfo$EdmondsonGrade==3)]="darkorange3"
  col.samples[which(sampleinfo$TissueType=="Tumor" & sampleinfo$EdmondsonGrade==4)]="red"
  
  names(col.samples)=sampleinfo$BiopsyID


  liver.samples=sampleinfo$BiopsyID[which(sampleinfo$TissueType=="Liver")]
  tumor.grade1=sampleinfo$BiopsyID[which(sampleinfo$TissueType=="Tumor" & sampleinfo$EdmondsonGrade==1)]
  tumor.grade2=sampleinfo$BiopsyID[which(sampleinfo$TissueType=="Tumor" & sampleinfo$EdmondsonGrade==2)]
  tumor.grade3=sampleinfo$BiopsyID[which(sampleinfo$TissueType=="Tumor" & sampleinfo$EdmondsonGrade==3)]
  tumor.grade4=sampleinfo$BiopsyID[which(sampleinfo$TissueType=="Tumor" & sampleinfo$EdmondsonGrade==4)]
  
########################################################################

  tpm=read.table(paste(pathExpression, annot, "/AllSamples_KallistoRawTPM.txt",sep=""), h=T, stringsAsFactors=F)

  tpm=tpm[,sampleinfo$BiopsyID]
  
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

pdf(file=paste(pathResults, "Figure1.pdf",sep=""),width=7, height=5.75)

m=matrix(c(rep(NA,16*8)), nrow=16)

for(i in 1:10){
  m[i,]=c(rep(1, 4), rep(2, 4))
}


for(i in 11:16){
  m[i,]=c(rep(3, 2), rep(4, 2), rep(5, 2), rep(6,2))
}

layout(m)

####################################################################

par(mar=c(6.5,4.5,2.1,1.1))

## protein-coding.genes

ylim=range(pca.pc$li[,2])
ylim[2]=ylim[2]+diff(ylim)/10

xlim=range(pca.pc$li[,1])
xlim[1]=xlim[1]-diff(xlim)/10

plot(pca.pc$li[,1], pca.pc$li[,2], pch=21, col="gray40", bg=col.samples[rownames(pca.pc$li)], xlab="", ylab="", axes=F, cex=1.25, ylim=ylim, xlim=xlim)
axis(side=1, mgp=c(3, 0.75, 0), cex.axis=1.15)
axis(side=2, mgp=c(3, 0.5, 0), cex.axis=1.15)
box()
mtext("protein-coding genes", side=3, line=0.35, cex=0.8)
mtext(paste("PC1 (", explained.pc[1],"% variance)",sep=""), side=1, line=2.45, cex=0.8)
mtext(paste("PC2 (", explained.pc[2],"% variance)",sep=""), side=2, line=2.35, cex=0.8)

legend("topleft", legend=c("liver", "tumor grade1", "tumor grade 2", "tumor grade 3", "tumor grade 4"), pt.bg=c("gray40", "gold", "orange", "darkorange3", "red"), bty="n", col="gray40", pch=21, inset=0.01)

mtext("A", side=3, at=-205, line=0.5, font=2, cex=1)

########################################################################

## lncRNAs


ylim=range(pca.lnc$li[,2])
ylim[2]=ylim[2]+diff(ylim)/10

xlim=range(pca.lnc$li[,1])
xlim[1]=xlim[1]-diff(xlim)/10

plot(pca.lnc$li[,1], pca.lnc$li[,2], pch=21, col="gray40", bg=col.samples[rownames(pca.lnc$li)], ylab="", xlab="",axes=F, cex=1.25, ylim=ylim, xlim=xlim)
axis(side=1, mgp=c(3, 0.75, 0), cex.axis=1.15)
axis(side=2, mgp=c(3, 0.5, 0), cex.axis=1.15)
box()
mtext("lncRNAs", side=3, line=0.35, cex=0.8)
mtext(paste("PC1 (", explained.lnc[1],"% variance)",sep=""), side=1, line=2.45, cex=0.8)
mtext(paste("PC2 (", explained.lnc[2],"% variance)",sep=""), side=2, line=2.35, cex=0.8)

mtext("B", side=3, at=-55, line=0.5, font=2, cex=1)

########################################################################

## protein-coding pc1

par(mar=c(2.2,4.5,0.5,0.5))

boxplot(pca.pc$li[liver.samples,1], pca.pc$li[tumor.grade1,1], pca.pc$li[tumor.grade2,1], pca.pc$li[tumor.grade3,1], pca.pc$li[tumor.grade4,1], horizontal=F, col="white", border=c("gray40", "gold", "orange", "darkorange3", "red"), outline=F, ylim=range(pca.pc$li[,1]))

mtext("PC1", side=2, line=2.45, cex=0.8)

mtext("sample types", side=1, line=1, at=7.5, cex=0.8)

mtext("C", side=3, at=-1.75, line=0.5, font=2, cex=1)

## protein-coding pc2

par(mar=c(2.2,4.15,0.5,0.4))

boxplot(pca.pc$li[liver.samples,2], pca.pc$li[tumor.grade1,2], pca.pc$li[tumor.grade2,2], pca.pc$li[tumor.grade3,2], pca.pc$li[tumor.grade4,2], horizontal=F, col="white", border=c("gray40", "gold", "orange", "darkorange3", "red"), outline=F, ylim=range(pca.pc$li[,2]))

mtext("PC2", side=2, line=2.45, cex=0.8)

mtext("D", side=3, at=-1.5, line=0.5, font=2, cex=1)

########################################################################

## lncRNA pc1

par(mar=c(2.2,4.5,0.5,0.5))

boxplot(pca.lnc$li[liver.samples,1], pca.lnc$li[tumor.grade1,1], pca.lnc$li[tumor.grade2,1], pca.lnc$li[tumor.grade3,1], pca.lnc$li[tumor.grade4,1], horizontal=F, col="white", border=c("gray40", "gold", "orange", "darkorange3", "red"), outline=F, ylim=range(pca.lnc$li[,1]))

mtext("PC1", side=2, line=2.45, cex=0.8)

mtext("sample types", side=1, line=1, at=7.5, cex=0.8)

mtext("E", side=3, at=-1.75, line=0.5, font=2, cex=1)

## protein-coding pc2

par(mar=c(2.2,4.15,0.5,0.4))

boxplot(pca.lnc$li[liver.samples,2], pca.lnc$li[tumor.grade1,2], pca.lnc$li[tumor.grade2,2], pca.lnc$li[tumor.grade3,2], pca.lnc$li[tumor.grade4,2], horizontal=F, col="white", border=c("gray40", "gold", "orange", "darkorange3", "red"), outline=F, ylim=range(pca.lnc$li[,2]))

mtext("PC2", side=2, line=2.45, cex=0.8)

mtext("F", side=3, at=-1.5, line=0.5, font=2, cex=1)


########################################################################

dev.off()

########################################################################
