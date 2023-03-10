####################################################################################

pathAnnot="../../data/ensembl_annotations/"
pathResults="../../results/differential_expression/"
pathPCHiC="../../results/PCHiC_interactions/"

release=109
annot="AllTranscripts_Ensembl109"
expdata="AllTranscripts_Ensembl109_noMT_norRNA_nohaplo"

maxFDR=1e-3
minFC=1.5

####################################################################################

## gene types

geneinfo=read.table(paste(pathAnnot, "GeneInfo_Ensembl",release,".txt",sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
rownames(geneinfo)=geneinfo$Gene.stable.ID

pc=geneinfo$Gene.stable.ID[which(geneinfo$Gene.type=="protein_coding")]
lnc=geneinfo$Gene.stable.ID[which(geneinfo$Gene.type=="lncRNA")]

####################################################################################

## contacts

contacts=read.table(paste(pathPCHiC, "PromoterPromoterInteractions_AllTranscripts_Ensembl109_WindowSize1000_unbaited_genesummary.txt", sep=""), h=T, stringsAsFactors=F)
contacts=contacts[which(contacts$Gene1%in%c(pc,lnc) & contacts$Gene2%in%c(pc,lnc) & contacts$MinDistance>50000 & contacts$MinDistance<2.5e6),]

####################################################################################

compute.prop.contacts <- function(deres, contacts, up.pc, up.lnc, down.pc, down.lnc){
    signif.genes=c(up.pc, up.lnc, down.pc, down.lnc)
    up.genes=c(up.pc, up.lnc)
    down.genes=c(down.pc, down.lnc)

    signif.contacts=lapply(signif.genes, function(x) unique(c(contacts$Gene2[which(contacts$Gene1==x)]))) ## look only at the genes contacted by the baited gene
    names(signif.contacts)=signif.genes

    prop.up.pc.has.up.contact=length(which(unlist(lapply(signif.contacts[up.pc], function(x) any(x%in%up.genes)))))/length(up.pc)
    prop.up.lnc.has.up.contact=length(which(unlist(lapply(signif.contacts[up.lnc], function(x) any(x%in%up.genes)))))/length(up.lnc)

    prop.down.pc.has.down.contact=length(which(unlist(lapply(signif.contacts[down.pc], function(x) any(x%in%down.genes)))))/length(down.pc)
    prop.down.lnc.has.down.contact=length(which(unlist(lapply(signif.contacts[down.lnc], function(x) any(x%in%down.genes)))))/length(down.lnc)

    results=list("pc.upup"=prop.up.pc.has.up.contact, "pc.downdown"=prop.down.pc.has.down.contact, "lnc.upup"=prop.up.lnc.has.up.contact, "lnc.downdown"=prop.down.lnc.has.down.contact)
    return(results)
}

####################################################################################

for(test in c("Tumor_vs_Liver", "EdmondsonGrade_34_vs_12")){
    deres=read.table(paste(pathResults,expdata, "/DifferentialExpression_",test,".txt", sep=""), h=T, stringsAsFactors=F, sep="\t")

    ## take only genes that are baited
    deres=deres[which(rownames(deres)%in%contacts$Gene1),]

    up.genes=rownames(deres)[which(deres$padj<maxFDR & deres$log2FoldChange>=log2(minFC))]
    down.genes=rownames(deres)[which(deres$padj<maxFDR & deres$log2FoldChange<=log2(1/minFC))]

    up.pc=intersect(up.genes, pc)
    down.pc=intersect(down.genes, pc)
    up.lnc=intersect(up.genes, lnc)
    down.lnc=intersect(down.genes, lnc)

    results=compute.prop.contacts(deres, contacts, up.pc, up.lnc, down.pc, down.lnc)

    all.results=list()
    all.results[["real"]]=unlist(results)

    ## we randomize significant genes

    for(i in 1:100){
        print(i)
        random.up.pc=sample(intersect(pc, rownames(deres)),  size=length(up.pc), re=FALSE)
        random.down.pc=sample(setdiff(intersect(pc, rownames(deres)), random.up.pc),  size=length(down.pc), re=FALSE)

        random.up.lnc=sample(intersect(lnc, rownames(deres)),  size=length(up.lnc), re=FALSE)
        random.down.lnc=sample(setdiff(intersect(lnc, rownames(deres)), random.up.lnc),  size=length(down.lnc), re=FALSE)

        random.results=compute.prop.contacts(deres, contacts, random.up.pc, random.up.lnc, random.down.pc, random.down.lnc)
        all.results[[paste("rand",i,sep="")]]=unlist(random.results)
    }

    all.results=t(as.data.frame(all.results))

    all.results=as.data.frame(all.results)

    write.table(all.results, file=paste(pathResults, expdata, "/ProportionDE_ContactedGenes_minDistance50kb_maxDistance2.5Mb_",test,"_maxFDR", maxFDR, "_minFC",minFC,".txt",sep=""), row.names=T, col.names=T, sep="\t")

}

####################################################################################

