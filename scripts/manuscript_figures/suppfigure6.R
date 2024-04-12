#############################################################################

if(!("pathFigures"%in%objects())){
    pathRData="../../data_for_publication/RData/"
    pathFigures="../../data_for_publication/main_figures/"

    library(vioplot)

    load=TRUE
    prepare=TRUE
}

#############################################################################

if(load){
    ## cited lncRNAs
    load(paste(pathRData, "data.PubMed.analysis.RData",sep=""))

    ## gene info
    load(paste(pathRData, "data.gene.info.RData", sep=""))

    ## differential expression
    load(paste(pathRData, "data.diffexp.RData",sep=""))

    maxFDR=0.05
    minFC=1.5

    load=FALSE
}

#############################################################################

if(prepare){

    lnc.cited.once=names(nb.citations.lnc)[which(nb.citations.lnc==1)]
    lnc.cited.more=names(nb.citations.lnc)[which(nb.citations.lnc>1)]
    lnc.cited.all=names(nb.citations.lnc)
    other.lnc=setdiff(lnc, lnc.cited.all)

    pc.cited.once=names(nb.citations.pc)[which(nb.citations.pc==1)]
    pc.cited.more=names(nb.citations.pc)[which(nb.citations.pc>1)]
    pc.cited.all=names(nb.citations.pc)
    other.pc=setdiff(pc, pc.cited.all)

    ## extract significant genes for the two DE analyses

    strong.upregulated=list()
    strong.downregulated=list()

    soft.upregulated=list()
    soft.downregulated=list()
    
    for(type in c("tnt.meric", "tnt.tcga")){
        this.de=get(paste("diffexp", type, sep="."))

        this.soft.up=rownames(this.de)[which(this.de$padj<maxFDR & this.de$log2FoldChange>0)]
        this.soft.down=rownames(this.de)[which(this.de$padj<maxFDR & this.de$log2FoldChange<0)]

        this.strong.up=rownames(this.de)[which(this.de$padj<maxFDR & this.de$log2FoldChange>log2(minFC))]
        this.strong.down=rownames(this.de)[which(this.de$padj<maxFDR & this.de$log2FoldChange<log2(1/minFC))]

        strong.upregulated[[type]]=this.strong.up
        strong.downregulated[[type]]=this.strong.down
        
        soft.upregulated[[type]]=this.soft.up
        soft.downregulated[[type]]=this.soft.down
    }

    highly.consistent.up=intersect(strong.upregulated[["tnt.meric"]], strong.upregulated[["tnt.tcga"]])
    highly.consistent.down=intersect(strong.downregulated[["tnt.meric"]], strong.downregulated[["tnt.tcga"]])
    highly.consistent=c(highly.consistent.up, highly.consistent.down)
    
    weakly.consistent.up=c(intersect(strong.upregulated[["tnt.meric"]], soft.upregulated[["tnt.tcga"]]), intersect(soft.upregulated[["tnt.meric"]], strong.upregulated[["tnt.tcga"]]))
    weakly.consistent.down=c(intersect(strong.downregulated[["tnt.meric"]], soft.downregulated[["tnt.tcga"]]), intersect(soft.downregulated[["tnt.meric"]], strong.downregulated[["tnt.tcga"]]))
    weakly.consistent=c(weakly.consistent.up, weakly.consistent.down)
    
    highly.inconsistent=c(intersect(strong.upregulated[["tnt.meric"]], strong.downregulated[["tnt.tcga"]]), intersect(strong.downregulated[["tnt.meric"]], strong.upregulated[["tnt.tcga"]]))
    weakly.inconsistent=c(intersect(strong.upregulated[["tnt.meric"]], soft.downregulated[["tnt.tcga"]]), intersect(strong.downregulated[["tnt.meric"]], soft.upregulated[["tnt.tcga"]]), intersect(soft.upregulated[["tnt.meric"]], strong.downregulated[["tnt.tcga"]]), intersect(soft.downregulated[["tnt.meric"]], strong.upregulated[["tnt.tcga"]]))
   
    
    
    prepare=FALSE
}

#############################################################################

pdf(file=paste(pathFigures, "SupplementaryFigure6.pdf", sep=""), width=6.85, height=7)

m=matrix(rep(NA, 2*10), nrow=2)

m[1,]=c(rep(1,7), rep(2,3))
m[1,]=c(rep(3,7), rep(4,3))

layout(m)

#############################################################################

for(gene

#############################################################################

dev.off()

#############################################################################
