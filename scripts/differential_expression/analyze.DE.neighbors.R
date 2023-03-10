####################################################################################

pathAnnot="../../data/ensembl_annotations/"
pathResults="../../results/differential_expression/"
pathGeneOverlaps="../../results/gene_overlaps/"

release=109
annot="AllTranscripts_Ensembl109"
expdata="AllTranscripts_Ensembl109_noMT_norRNA_nohaplo"

maxFDR=1e-5
minFC=1.5

####################################################################################

## gene types

geneinfo=read.table(paste(pathAnnot, "GeneInfo_Ensembl",release,".txt",sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
rownames(geneinfo)=geneinfo$Gene.stable.ID

pc=geneinfo$Gene.stable.ID[which(geneinfo$Gene.type=="protein_coding")]
lnc=geneinfo$Gene.stable.ID[which(geneinfo$Gene.type=="lncRNA")]

####################################################################################

## overlapping genes

geneov=read.table(paste(pathGeneOverlaps, "OverlappingGenes_", annot, ".txt",sep=""), h=T, stringsAsFactors=F)
geneov=geneov[which(geneov$GeneID%in%c(pc,lnc) & geneov$NeighborID%in%c(pc,lnc)),]

####################################################################################

## neighboring genes

neighbors=read.table(paste(pathGeneOverlaps, "NeighboringGenes_50kb_",annot,".txt",sep=""), h=T, stringsAsFactors=F)
neighbors=neighbors[which(neighbors$GeneID%in%c(pc,lnc) & neighbors$NeighborID%in%c(pc,lnc)),]

####################################################################################

compute.prop.neighbors <- function(deres, neighbors, up.pc, up.lnc, down.pc, down.lnc){
    signif.genes=c(up.pc, up.lnc, down.pc, down.lnc)
    up.genes=c(up.pc, up.lnc)
    down.genes=c(down.pc, down.lnc)

    signif.neighbors=lapply(signif.genes, function(x) neighbors$NeighborID[which(neighbors$GeneID==x)])
    names(signif.neighbors)=signif.genes

    prop.up.pc.has.up.neighbor=length(which(unlist(lapply(signif.neighbors[up.pc], function(x) any(x%in%up.genes)))))/length(up.pc)
    prop.up.lnc.has.up.neighbor=length(which(unlist(lapply(signif.neighbors[up.lnc], function(x) any(x%in%up.genes)))))/length(up.lnc)

    prop.down.pc.has.down.neighbor=length(which(unlist(lapply(signif.neighbors[down.pc], function(x) any(x%in%down.genes)))))/length(down.pc)
    prop.down.lnc.has.down.neighbor=length(which(unlist(lapply(signif.neighbors[down.lnc], function(x) any(x%in%down.genes)))))/length(down.lnc)

    results=list("pc.upup"=prop.up.pc.has.up.neighbor, "pc.downdown"=prop.down.pc.has.down.neighbor, "lnc.upup"=prop.up.lnc.has.up.neighbor, "lnc.downdown"=prop.down.lnc.has.down.neighbor)
    return(results)
}

####################################################################################

for(test in c("Tumor_vs_Liver", "EdmondsonGrade_34_vs_12")){
    deres=read.table(paste(pathResults,expdata, "/DifferentialExpression_",test,".txt", sep=""), h=T, stringsAsFactors=F, sep="\t")

    up.genes=rownames(deres)[which(deres$padj<maxFDR & deres$log2FoldChange>=log2(minFC))]
    down.genes=rownames(deres)[which(deres$padj<maxFDR & deres$log2FoldChange<=log2(1/minFC))]

    up.pc=intersect(up.genes, pc)
    down.pc=intersect(down.genes, pc)
    up.lnc=intersect(up.genes, lnc)
    down.lnc=intersect(down.genes, lnc)

    results=compute.prop.neighbors(deres, neighbors, up.pc, up.lnc, down.pc, down.lnc)

    all.results=list()
    all.results[["real"]]=unlist(results)

    ## we randomize significant genes

    for(i in 1:100){
        print(i)
        random.up.pc=sample(intersect(pc, rownames(deres)),  size=length(up.pc), re=FALSE)
        random.down.pc=sample(setdiff(intersect(pc, rownames(deres)), random.up.pc),  size=length(down.pc), re=FALSE)

        random.up.lnc=sample(intersect(lnc, rownames(deres)),  size=length(up.lnc), re=FALSE)
        random.down.lnc=sample(setdiff(intersect(lnc, rownames(deres)), random.up.lnc),  size=length(down.lnc), re=FALSE)

        random.results=compute.prop.neighbors(deres, neighbors, random.up.pc, random.up.lnc, random.down.pc, random.down.lnc)
        all.results[[paste("rand",i,sep="")]]=unlist(random.results)
    }

    all.results=as.data.frame(all.results)

}

####################################################################################

