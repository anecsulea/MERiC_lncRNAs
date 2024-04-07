###########################################################################

pathDocs="../../docs/"
pathSNPAnalysis="../../results/SNP_analysis/"
pathRData="../../data_for_publication/RData/"

###########################################################################

load(paste(pathRData, "data.sample.info.MERiC.RData",sep=""))

###########################################################################

## check percentage of shared SNPs

samples.per.patient=lapply(unique(nontumor.samples$Patient_ID), function(x) c(nontumor.samples$biopsyID[which(nontumor.samples$Patient_ID==x)], tumor.samples$tumor_biopsyID[which(tumor.samples$Patient_ID==x)]))
names(samples.per.patient)=unique(nontumor.samples$Patient_ID)

shared.alleles=read.table("../../results/SNP_analysis/ProportionSharedAlleles_BiallelicSNPs_MinCoverage50_NoRepeats_GATK.txt", h=F, stringsAsFactors=F)

pc.shared.alleles=list()

for(patient in names(samples.per.patient)){
    samples=samples.per.patient[[patient]]

    pc.shared.alleles[[patient]]=c()

    for(i in 1:(length(samples)-1)){
        sample1=samples[i]
        for(j in (i+1):length(samples)){
            sample2=samples[j]
            pc.shared=shared.alleles[which(shared.alleles$V1==sample1 & shared.alleles$V2==sample2)[1],3]

            pc.shared.alleles[[patient]]=c(pc.shared.alleles[[patient]], pc.shared)
        }
    }
}

###########################################################################
