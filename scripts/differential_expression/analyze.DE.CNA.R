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

####################################################################################

tpm=read.table(paste(pathExpression, expdata, "/AllSamples_KallistoRawTPM.txt",sep=""), h=T, stringsAsFactors=F, sep="\t")

####################################################################################

## sample types
sampleinfo=read.table(paste(pathSampleInfo, "AnalyzedSamples.txt",sep=""), h=T,stringsAsFactors=F, sep="\t", quote="\"")
liver.samples=sampleinfo$BiopsyID[which(sampleinfo$TissueType=="Liver")]
tumor.samples=sampleinfo$BiopsyID[which(sampleinfo$TissueType=="Tumor")]

####################################################################################
## average expression

avg.exp.liver=apply(tpm[,liver.samples],1,mean)
avg.exp.tumor=apply(tpm[,tumor.samples],1,mean)

####################################################################################

## gene types

geneinfo=read.table(paste(pathAnnot, "GeneInfo_Ensembl",release,".txt",sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
rownames(geneinfo)=geneinfo$Gene.stable.ID

geneinfo=geneinfo[which(geneinfo$Chromosome.scaffold.name%in%c(as.character(1:22), "X", "Y")),]

pc=geneinfo$Gene.stable.ID[which(geneinfo$Gene.type=="protein_coding")]
lnc=geneinfo$Gene.stable.ID[which(geneinfo$Gene.type=="lncRNA")]
pseudo=geneinfo$Gene.stable.ID[which(geneinfo$Gene.type%in%c("transcribed_unitary_pseudogene", "transcribed_unprocessed_pseudogene"))]

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

compute.prop.cna <- function(recurrently.amp, recurrently.del, up.pc, up.lnc, down.pc, down.lnc){
    signif.genes=c(up.pc, up.lnc, down.pc, down.lnc)
    up.genes=c(up.pc, up.lnc)
    down.genes=c(down.pc, down.lnc)

    prop.up.pc.amp=length(intersect(up.pc, recurrently.amp))/length(up.pc)
    prop.up.lnc.amp=length(intersect(up.lnc, recurrently.amp))/length(up.lnc)

    prop.down.pc.del=length(intersect(down.pc, recurrently.del))/length(down.pc)
    prop.down.lnc.del=length(intersect(down.lnc, recurrently.del))/length(down.lnc)

    prop.up.pc.del=length(intersect(up.pc, recurrently.del))/length(up.pc)
    prop.up.lnc.del=length(intersect(up.lnc, recurrently.del))/length(up.lnc)

    prop.down.pc.amp=length(intersect(down.pc, recurrently.amp))/length(down.pc)
    prop.down.lnc.amp=length(intersect(down.lnc, recurrently.amp))/length(down.lnc)

    results=list("pc.up.amp"=prop.up.pc.amp, "pc.up.del"=prop.up.pc.del,  "lnc.up.amp"=prop.up.lnc.amp, "lnc.up.del"=prop.up.lnc.del, "pc.down.amp"=prop.down.pc.amp, "pc.down.del"=prop.down.pc.del,  "lnc.down.amp"=prop.down.lnc.amp, "lnc.down.del"=prop.down.lnc.del)
    
    return(results)
}

####################################################################################

for(test in c("Tumor_vs_Liver", "EdmondsonGrade_34_vs_12")){
  
  deres=read.table(paste(pathResults,expdata, "/DifferentialExpression_",test,".txt", sep=""), h=T, stringsAsFactors=F, sep="\t")
  
  up.genes=rownames(deres)[which(deres$padj<maxFDR & deres$log2FoldChange>=log2(minFC))]
  down.genes=rownames(deres)[which(deres$padj<maxFDR & deres$log2FoldChange<=log2(1/minFC))]
  
  up.pc=intersect(pc, up.genes)
  down.pc=intersect(pc, down.genes)
  tested.pc=intersect(pc, rownames(deres)[which(!is.na(deres$padj))])
  
  up.lnc=intersect(lnc, up.genes)
  down.lnc=intersect(lnc, down.genes)
  tested.lnc=intersect(lnc, rownames(deres)[which(!is.na(deres$padj))])

  results=compute.prop.cna(recurrently.amp, recurrently.del, up.pc, up.lnc, down.pc, down.lnc)
  
  all.results=list()
  all.results[["real"]]=unlist(results)
  
  ## we randomize significant genes
  
  for(i in 1:100){
    print(i)
    random.up.pc=sample(intersect(pc, rownames(deres)),  size=length(up.pc), re=FALSE)
    random.down.pc=sample(setdiff(intersect(pc, rownames(deres)), random.up.pc),  size=length(down.pc), re=FALSE)
    
    random.up.lnc=sample(intersect(lnc, rownames(deres)),  size=length(up.lnc), re=FALSE)
    random.down.lnc=sample(setdiff(intersect(lnc, rownames(deres)), random.up.lnc),  size=length(down.lnc), re=FALSE)
    
    random.results=compute.prop.cna(recurrently.amp, recurrently.del, random.up.pc, random.up.lnc, random.down.pc, random.down.lnc)
    all.results[[paste("rand",i,sep="")]]=unlist(random.results)
  }
  
  all.results=t(as.data.frame(all.results))
  
  all.results=as.data.frame(all.results)

  write.table(all.results, file=paste(pathResults, expdata, "/ProportionDE_CNA_",test,"_maxFDR", maxFDR, "_minFC",minFC,".txt",sep=""), row.names=T, col.names=T, sep="\t")
}


####################################################################################
