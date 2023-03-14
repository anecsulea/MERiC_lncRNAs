##############################################################################

pathExpression="../../results/expression_estimation/"
pathDifferentialExpression="../../results/differential_expression/"
pathAnnot="../../data/ensembl_annotations/"
pathSampleInfo="../../results/sample_info/"
pathResults="../../results/figures/main_figures/"

annot="AllTranscripts_Ensembl109_noMT_norRNA_nohaplo"

set.seed(19)

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

  
}

########################################################################
#######################################################################

pdf(file=paste(pathResults, "PCA_gene_types.pdf",sep=""),width=9,height=5)

p
par(mar=c(4.1,2.1,1.1,1.1))

plot(pca.pc$li[,1], pca.pc$li[,2], pch=20, col=col.samples[rownames(pca.pc$li)])

plot(pca.lnc$li[,1], pca.lnc$li[,2], pch=20, col=col.samples[rownames(pca.lnc$li)])

plot(pca.pseudo$li[,1], pca.pseudo$li[,2], pch=20, col=col.samples[rownames(pca.pseudo$li)])
dev.off()

########################################################################
