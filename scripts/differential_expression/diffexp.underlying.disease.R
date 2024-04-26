########################################################################

library(DESeq2)
library(BiocParallel)
library(tximport)

register(MulticoreParam(4))

########################################################################

pathExpression="../../results/expression_estimation/"
pathDifferentialExpression="../../results/differential_expression/"
pathRData="../../data_for_publication/RData/"

annot="AllTranscripts_Ensembl109_noMT_norRNA_nohaplo"

########################################################################

load(paste(pathRData, "data.sample.info.MERiC.RData",sep=""))

rownames(tumor.samples)=tumor.samples$tumor_biopsyID

########################################################################

## get read counts for each gene

read.counts=read.table(paste(pathExpression, annot, "/AllSamples_KallistoEstimatedCounts.txt",sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")

total.estimated.counts=apply(read.counts,2,sum)

tumor.samples$total.estimated.counts=total.estimated.counts[tumor.samples$tumor_biopsyID]

## order samples by read count and Edmondson degree
tumor.samples=tumor.samples[order(tumor.samples$total.estimated.counts, decreasing=T),]

tumor.samples=tumor.samples[order(tumor.samples$edmondson, decreasing=T),]

## select one sample per patient

dupli=which(duplicated(tumor.samples$Patient_ID))

if(length(dupli)>0){
    tumor.samples=tumor.samples[-dupli,]
}

########################################################################

## get tx2gene info

samples=tumor.samples$tumor_biopsyID

sinfo=read.table(paste(pathExpression, annot, "/",samples[1],"/abundance.tsv",sep=""), h=T, stringsAsFactors=F)
geneid=unlist(lapply(sinfo$target_id, function(x) unlist(strsplit(x, split=":"))[1]))

tx2gene=data.frame("txid"=sinfo$target_id, "geneid"=geneid)

########################################################################

disease.synonyms=c("ALD", "NAFLD", "HepB","HepC")
names(disease.synonyms)= c("Alcoholic liver disease", "Non-alcoholic liver disease", "Hepatitis B", "Hepatitis C")

for(disease in c("Alcoholic liver disease", "Non-alcoholic liver disease", "Hepatitis B", "Hepatitis C")){

    diseased.samples=tumor.samples$tumor_biopsyID[grep(disease, tumor.samples$underlying_liver_disease)]
    no.disease=tumor.samples$tumor_biopsyID[which(tumor.samples$underlying_liver_disease=="No liver disease")]
    other.samples=setdiff(tumor.samples$tumor_biopsyID, diseased.samples)

    for(control in c("no.disease", "other.samples")){

        control.samples=get(control)

        all.samples=c(diseased.samples, control.samples)

        print(paste(disease, control))
        print(paste(length(diseased.samples), "disease samples"))
        print(paste(length(control.samples), "control.samples"))
                
        ## assemble kallisto counts
        
        files=paste(pathExpression, "/", annot,"/", all.samples, "/abundance.h5", sep="")
        names(files)=all.samples
        
        txi.kallisto <- tximport(files, type = "kallisto", tx2gene = tx2gene)
        
        ## prepare design

        disease.factor=c(rep("disease", length(diseased.samples)), rep("control", length(control.samples)))

        this.edmondson=tumor.samples[all.samples, "edmondson"]
        
        eg=rep(NA, length(all.samples))
        
        eg[which(this.edmondson%in%c(1,2))]="12"
        eg[which(this.edmondson%in%c(3,4))]="34"
        
        colData=data.frame("Sex"=as.factor(tumor.samples[all.samples, "sex"]), "EdmondsonGrade"=as.factor(eg), "Disease"=as.factor(disease.factor))

        print(colData)
        
        ## run DESeq2
        
        dds=DESeqDataSetFromTximport(txi.kallisto, colData=colData, design = ~Sex+EdmondsonGrade+Disease)

        dds=DESeq(dds, test="Wald",  minReplicatesForReplace=10, parallel=T)

        ## get test results

        results=results(dds, contrast=c("Disease", "disease", "control"))

        results=lfcShrink(dds, coef="Disease_disease_vs_control", res=results, type="apeglm")
        
        results=as.data.frame(results)
        
        results=results[order(results$padj),]
        
        ## write results
        
        syn=disease.synonyms[disease]
        
        write.table(results, file=paste(pathDifferentialExpression, annot, "/DifferentialExpression_", syn,"_vs_",control,".txt",sep=""), row.names=T, col.names=T, sep="\t", quote=F)

    }
}

########################################################################
