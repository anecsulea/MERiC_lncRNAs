###############################################################################

annot="AllTranscripts_Ensembl109_noMT_norRNA_nohaplo"

pathExpression=paste("../../results/expression_estimation/", annot,sep="")
pathResults=paste("../../results/co_expression_networks/",annot,sep="")
pathEnsembl="../../data/ensembl_annotations/"

release=109

options(stringsAsFactors=F)

###############################################################################

exp=read.table(paste(pathExpression, "/AllSamples_KallistoRawTPM.txt",sep=""), h=T, stringsAsFactors=F, sep="\t")

###############################################################################

info=read.table(paste(pathEnsembl, "GeneInfo_Ensembl",release,".txt",sep=""), h=T, stringsAsFactors=F,sep="\t", quote="\"")
info=info[which(info$Chromosome.scaffold.name%in%c(as.character(1:22), "X", "Y")),] ## well-assembled chromosomes
rownames(info)=info$Gene.stable.ID

####################################################################

pc=info$Gene.stable.ID[which(info$Gene.type=="protein_coding")]
lnc=info$Gene.stable.ID[which(info$Gene.type=="lncRNA")]

all.genes=intersect(c(pc, lnc), rownames(exp))

####################################################################

exp=exp[all.genes,]

maxexp=apply(exp,1, max)

selected.genes=all.genes[which(maxexp>1)]

####################################################################

aracne.exp=data.frame("gene"=selected.genes)
aracne.exp=cbind(aracne.exp, exp[selected.genes,])

write.table(aracne.exp, file=paste(pathResults, "/expression_data.txt",sep=""), row.names=F, col.names=T, sep="\t", quote=F)
write.table(selected.genes, file=paste(pathResults, "/all_genes.txt",sep=""), row.names=F, col.names=T, sep="\t", quote=F)
write.table(intersect(selected.genes, lnc), file=paste(pathResults, "/lncRNAs.txt",sep=""), row.names=F, col.names=T, sep="\t", quote=F)

####################################################################
