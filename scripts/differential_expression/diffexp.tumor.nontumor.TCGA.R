########################################################################

set.seed(19)

options(stringsAsFactors=F)

########################################################################

library(DESeq2)
library(BiocParallel)
library(tximport)

register(MulticoreParam(8))

########################################################################

pathExpression="../../results/expression_estimation_TCGA/"
pathDifferentialExpression="../../results/differential_expression_TCGA/"
pathTCGA="../../data/TCGA/"
pathRData="../../data_for_publication/RData/"

annot="AllTranscripts_Ensembl109_noMT_norRNA_nohaplo"

########################################################################
## read sample info

sampleinfo=read.table(paste(pathTCGA, "all_info_paired_samples.txt", sep=""), h=T, stringsAsFactors=F, sep="\t")
sampleinfo$case_id=as.character(sampleinfo$case_id)
rownames(sampleinfo)=sampleinfo$id

sampleinfo$sample_type[which(sampleinfo$sample_type=="Solid Tissue Normal")]="Non-Tumor"
sampleinfo$sample_type[which(sampleinfo$sample_type=="Primary Tumor")]="Tumor"

########################################################################

## get tx2gene info

samples=sampleinfo$id

sinfo=read.table(paste(pathExpression, annot, "/",samples[1],"/abundance.tsv",sep=""), h=T, stringsAsFactors=F)
geneid=unlist(lapply(sinfo$target_id, function(x) unlist(strsplit(x, split=":"))[1]))

tx2gene=data.frame("txid"=sinfo$target_id, "geneid"=geneid)

########################################################################

## assemble kallisto counts

files=paste(pathExpression, "/", annot,"/", samples, "/abundance.h5", sep="")
names(files)=samples

txi.kallisto <- tximport(files, type = "kallisto", tx2gene = tx2gene)

########################################################################

## tissue and patient factor, all samples

selected.samples=sampleinfo$id
tissue.factor=as.factor(sampleinfo[, "sample_type"])
patient.factor=as.factor(as.character(sampleinfo[, "case_id"]))

colData=data.frame("tissue"=tissue.factor, "patient"=patient.factor)

## select samples in kallisto tximport table, same order

txi.kallisto$counts=txi.kallisto$counts[,selected.samples]
txi.kallisto$abundance=txi.kallisto$abundance[,selected.samples]
txi.kallisto$length=txi.kallisto$length[,selected.samples]

########################################################################

dds=DESeqDataSetFromTximport(txi.kallisto, colData=colData, design = ~tissue+patient)
dds=DESeq(dds, test="Wald",  minReplicatesForReplace=100, parallel=T)

## get DE genes

res=results(dds, name="tissue_Tumor_vs_Non.Tumor")
res.corrected=lfcShrink(dds, coef="tissue_Tumor_vs_Non.Tumor", res=res, type="apeglm")

res.corrected=res.corrected[order(res.corrected$pvalue),]
res.corrected=as.data.frame(res.corrected)

res.corrected$GeneID=rownames(res.corrected)
res.corrected=res.corrected[,c("GeneID", setdiff(colnames(res.corrected), "GeneID"))]

## write output table

write.table(res.corrected, file=paste(pathDifferentialExpression, annot, "/DifferentialExpressionResults_Tumor_vs_NonTumor.txt", sep=""), row.names=F, col.names=T, sep="\t", quote=F)

########################################################################
