##########################################################################

pathOverlaps="../../results/gene_overlaps/"
pathRData="../../data_for_publication/RData/"

annot="AllTranscripts_Ensembl109"

##########################################################################

## genic overlaps
gene.overlaps=read.table(paste(pathOverlaps, "OverlappingGenes_", annot,".txt",sep=""), h=T, stringsAsFactors=F)

## bidirectional promoters 1kb

biprom1kb=read.table(paste(pathOverlaps, "DistanceAntisenseTSS_MaxDist1kb_AllTranscripts_Ensembl109.txt", sep=""), h=T, stringsAsFactors=F)
biprom1kb=biprom1kb[, c("GeneID", "GenesCloseTSS")]

biprom5kb=read.table(paste(pathOverlaps, "DistanceAntisenseTSS_MaxDist5kb_AllTranscripts_Ensembl109.txt", sep=""), h=T, stringsAsFactors=F)
biprom5kb=biprom5kb[, c("GeneID", "GenesCloseTSS")]

neighbors50kb=read.table(paste(pathOverlaps, "NeighboringGenes_50kb_AllTranscripts_Ensembl109.txt", sep=""), h=T, stringsAsFactors=F)

##########################################################################

save(list=c("gene.overlaps", "biprom1kb", "biprom5kb", "neighbors50kb"), file=paste(pathRData, "data.gene.overlaps.RData",sep=""))

##########################################################################
