########################################################################

pathResults="../../results/splicing_analysis/"
pathDocs="../../docs/"

annot="AllTranscripts_Ensembl109"

########################################################################

tumor.samples=read.table(paste(pathDocs, "TumorSamples_Ng2022.txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
tumor.samples=tumor.samples[which(tumor.samples$tumor_biopsyID!=""),]

###########################################################################

liver.samples=read.table(paste(pathDocs, "LiverSamples_Ng2022.txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
liver.samples=liver.samples[which(liver.samples$biopsyID!="POOL"),]

###########################################################################

## select one sample per patient

dupli=which(duplicated(tumor.samples$Patient_ID))

if(length(dupli)>0){
    tumor.samples=tumor.samples[-dupli,]
}

###########################################################################

sample.info.t=data.frame("BiopsyID"=tumor.samples$tumor_biopsyID, "TissueType"=rep("Tumor", nrow(tumor.samples)), "Sex"=tumor.samples$sex, "AgeAtBiopsy"=tumor.samples$age_at_biopsy, "EdmondsonGrade"=tumor.samples$edmondson, "Cirrhosis"=tumor.samples$cirrhosis, "Diseases"=tumor.samples$underlying_liver_disease)

sample.info.l=data.frame("BiopsyID"=liver.samples$biopsyID, "TissueType"=rep("Liver", nrow(liver.samples)), "Sex"=liver.samples$sex, "AgeAtBiopsy"=liver.samples$age_at_biopsy, "EdmondsonGrade"=rep(NA, nrow(liver.samples)), "Cirrhosis"=rep(NA, nrow(liver.samples)), "Diseases"=rep(NA, nrow(liver.samples)))

sample.info=rbind(sample.info.t, sample.info.l)

samples=sample.info$BiopsyID

########################################################################

all.results=list()

for(i in 1:length(samples)){
  sample=samples[i]

  if(file.exists(paste(pathResults, sample, "/CoverageGenes_IntronBlocks_", annot, ".txt", sep="")) & file.exists(paste(pathResults, sample, "/CoverageGenes_ExonBlocks_", annot, ".txt", sep=""))){
      print(sample)

      introns=read.table(paste(pathResults, sample, "/CoverageGenes_IntronBlocks_", annot, ".txt", sep=""), h=T, stringsAsFactors=F, sep="\t")
      rownames(introns)=introns$GeneID

      exons=read.table(paste(pathResults, sample, "/CoverageGenes_ExonBlocks_", annot, ".txt", sep=""), h=T, stringsAsFactors=F, sep="\t")
      rownames(exons)=exons$GeneID

      common=intersect(exons$GeneID, introns$GeneID)

      ratio=introns[common, "MeanCoverageSense"]/(exons[common, "MeanCoverageSense"]+introns[common, "MeanCoverageSense"])

      names(ratio)=common

      results=data.frame("GeneID"=common, "CoverageExons"=exons[common, "MeanCoverageSense"], "CoverageIntrons"=introns[common, "MeanCoverageSense"], "Ratio"=ratio, stringsAsFactors=F)

      write.table(results, file=paste(pathResults, sample, "/IntronExonRatio_", annot, ".txt", sep=""), row.names=F, col.names=T, sep="\t", quote=F)

      if(length(all.results)==0){
          all.results[[sample]]=ratio
      } else{
          ## we reorder results
          firstsample=samples[1]
          ids=names(all.results[[firstsample]])
          all.results[[sample]]=ratio[ids]
      }
  }
}

ids=names(all.results[[samples[1]]])
all.results=as.data.frame(all.results)
all.results$GeneID=ids
all.results=all.results[,c("GeneID", setdiff(colnames(all.results), "GeneID"))]

write.table(all.results, file=paste(pathResults, "/AllSamples_IntronExonRatio_", annot, ".txt", sep=""), row.names=F, col.names=T, sep="\t", quote=F)

########################################################################

