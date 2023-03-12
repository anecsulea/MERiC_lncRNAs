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
minFC=2

####################################################################################

tpm=read.table(paste(pathExpression, expdata, "/AllSamples_KallistoRawTPM.txt",sep=""), h=T, stringsAsFactors=F, sep="\t")

####################################################################################

sampleinfo=read.table(paste(pathSampleInfo, "AnalzedSamples.txt",sep=""), h=T,stringsAsFactors=F, sep="\t", quote="\"")

####################################################################################

## gene types

geneinfo=read.table(paste(pathAnnot, "GeneInfo_Ensembl",release,".txt",sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
rownames(geneinfo)=geneinfo$Gene.stable.ID

geneinfo=geneinfo[which(geneinfo$Chromosome.scaffold.name%in%c(as.character(1:22), "X", "Y")),]

pc=geneinfo$Gene.stable.ID[which(geneinfo$Gene.type=="protein_coding")]
lnc=geneinfo$Gene.stable.ID[which(geneinfo$Gene.type=="lncRNA")]

####################################################################################

test="Tumor_vs_Liver"

deres=read.table(paste(pathResults,expdata, "/DifferentialExpression_",test,".txt", sep=""), h=T, stringsAsFactors=F, sep="\t")

up.genes=rownames(deres)[which(deres$padj<maxFDR & deres$log2FoldChange>=log2(minFC))]
down.genes=rownames(deres)[which(deres$padj<maxFDR & deres$log2FoldChange<=log2(1/minFC))]

####################################################################################

cna=read.table(paste(pathCNA, "GeneOverlap_",annot,"_CNA_AllExonv6_CRE.txt",sep=""), h=T, stringsAsFactors=F,sep="\t")
cna$biotype=geneinfo[cna$GeneID, "Gene.type"]
cna=cna[which(cna$Sample%in%sampleinfo$BiopsyID),]
cna=cna[which(cna$GeneID%in%rownames(deres)),]

amp=cna[which(cna$Type=="amp"),]
del=cna[which(cna$Type=="del"),]

nb.samples.amp=tapply(amp$Sample, as.factor(amp$GeneID), function(x) length(unique(x)))
names(nb.samples.amp)=levels(as.factor(amp$GeneID))

nb.samples.del=tapply(del$Sample, as.factor(del$GeneID), function(x) length(unique(x)))
names(nb.samples.del)=levels(as.factor(del$GeneID))

recurrently.amp=names(nb.samples.amp)[which(nb.samples.amp>=5)]
recurrently.del=names(nb.samples.del)[which(nb.samples.del>=5)]

####################################################################################
