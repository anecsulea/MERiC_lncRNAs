#########################################################################################

pathAnnot="../../data/ensembl_annotations/"
pathResults="../../results/differential_expression/"
pathExpression="../../results/expression_estimation/"
pathSampleInfo="../../results/sample_info/"
pathDocs="../../docs/"
pathCNA="../../results/CNA_analysis/"

release=109
annot="AllTranscripts_Ensembl109"
expdata="AllTranscripts_Ensembl109_noMT_norRNA_nohaplo"

maxFDR=1e-3
minFC=1.5

binsize=10e6

###########################################################################

geneinfo=read.table(paste(pathAnnot, "GeneInfo_Ensembl",release,".txt",sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
rownames(geneinfo)=geneinfo$Gene.stable.ID

geneinfo=geneinfo[which(geneinfo$Chromosome.scaffold.name%in%c(as.character(1:22), "X", "Y")),]

pc=geneinfo$Gene.stable.ID[which(geneinfo$Gene.type=="protein_coding")]
lnc=geneinfo$Gene.stable.ID[which(geneinfo$Gene.type=="lncRNA")]
pseudo=geneinfo$Gene.stable.ID[which(geneinfo$Gene.type%in%c("transcribed_unitary_pseudogene", "transcribed_unprocessed_pseudogene"))]

###########################################################################

chrsizes=read.table(paste(pathAnnot, "chr_sizes.txt",sep=""))

###########################################################################

colnames(geneinfo)[3]="Chr"
geneinfo$MidPos=(geneinfo$Gene.start..bp.+geneinfo$Gene.end..bp.)/2

####################################################################################

## sample types
sampleinfo=read.table(paste(pathSampleInfo, "AnalyzedSamples.txt",sep=""), h=T,stringsAsFactors=F, sep="\t", quote="\"")
liver.samples=sampleinfo$BiopsyID[which(sampleinfo$TissueType=="Liver")]
tumor.samples=sampleinfo$BiopsyID[which(sampleinfo$TissueType=="Tumor")]

###########################################################################

cna=read.table(paste(pathCNA, "GeneOverlap_",annot,"_CNA_AllExonv6_CRE.txt",sep=""), h=T, stringsAsFactors=F,sep="\t")
cna$biotype=geneinfo[cna$GeneID, "Gene.type"]
cna=cna[which(cna$Sample%in%sampleinfo$BiopsyID),]

amp=cna[which(cna$Type=="amp"),]
del=cna[which(cna$Type=="del"),]

nb.samples.amp=tapply(amp$Sample, as.factor(amp$GeneID), function(x) length(unique(x)))
names(nb.samples.amp)=levels(as.factor(amp$GeneID))

nb.samples.del=tapply(del$Sample, as.factor(del$GeneID), function(x) length(unique(x)))
names(nb.samples.del)=levels(as.factor(del$GeneID))

recurrently.amp=names(nb.samples.amp)[which(nb.samples.amp>=5)]
recurrently.del=names(nb.samples.del)[which(nb.samples.del>=5)]

####################################################################################

for(test in c("Tumor_vs_Liver", "EdmondsonGrade_34_vs_12")){
  
  deres=read.table(paste(pathResults,expdata, "/DifferentialExpression_",test,".txt", sep=""), h=T, stringsAsFactors=F, sep="\t")
  tested.genes=rownames(deres)[which(!is.na(deres$padj))]
  up.genes=rownames(deres)[which(deres$padj<maxFDR & deres$log2FoldChange>=log2(minFC))]
  down.genes=rownames(deres)[which(deres$padj<maxFDR & deres$log2FoldChange<=log2(1/minFC))]
  
  for(chr in c("17")){
    chrlen=chrsizes[which(chrsizes[,1]==chr),2]

    tested.thischr=tested.genes[which(geneinfo[tested.genes, "Chr"]==chr)]
    tested.thispos=geneinfo[tested.thischr, "MidPos"]
    names(tested.thispos)=tested.thischr

    bins=cut(tested.thispos, breaks=seq(from=0, to=chrlen, by=binsize))
    bin.midpos=seq(from=binsize, to=chrlen, by=binsize)-binsize/2
    bin.nbtested=tapply(tested.thischr, bins, function(x) length(x))
    bin.propup=tapply(tested.thischr, bins, function(x) length(which(x%in%up.genes))/length(x))
    bin.propdown=tapply(tested.thischr, bins, function(x) length(which(x%in%down.genes))/length(x))

    bin.propamp=tapply(tested.thischr, bins, function(x) length(which(x%in%recurrently.amp))/length(x))
    bin.propdel=tapply(tested.thischr, bins, function(x) length(which(x%in%recurrently.del))/length(x))

    par(mfrow=c(2,1))
    plot(bin.midpos, bin.propup)
    plot(bin.midpos, bin.propamp)
       
    plot(1,type="n", xlim=c(1, maxpos), ylim=c(-0.5,1.5), axes=F, xlab="", ylab="")
    abline(h=0.5)
    points(tested.thispos, rep(0.5, length(tested.thispos)))
    points(tested.thispos[intersect(tested.thischr, up.genes)], rep(1, length(intersect(tested.thischr, up.genes))), col="red", pch=20)
    points(tested.thispos[intersect(tested.thischr, down.genes)], rep(0, length(intersect(tested.thischr, down.genes))), col="blue", pch=20)
  }
}




###########################################################################

