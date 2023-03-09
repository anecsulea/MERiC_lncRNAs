########################################################################

pathExpression="../../results/expression_estimation/"

########################################################################

set.seed(19)

source("normalization.R")

########################################################################

library(tximport)

########################################################################

for(annot in c("AllTranscripts_Ensembl109_noMT_norRNA_nohaplo")){

    if(file.exists(paste(pathExpression,  annot, "/AllSamples_KallistoRawTPM.txt",sep=""))){
        print(paste("already done",annot))
    } else{

        samples=system(paste("ls ",pathExpression,  annot, " | grep -v txt",sep=""), intern=T)

        ## tx2gene table

        info=read.table(paste(pathExpression, annot, "/",samples[1],"/abundance.tsv", sep=""), h=T, stringsAsFactors=F, sep="\t")
        tx=info$target_id
        gene=unlist(lapply(tx, function(x) unlist(strsplit(x, split=":"))[1]))
        tx2gene=data.frame("tx"=tx, "gene"=gene, stringsAsFactors=F)

        ###############################################################################

        files=paste(pathExpression, annot, "/", samples, "/abundance.h5", sep="")
        names(files)=samples

        txi.kallisto <- tximport(files, type = "kallisto", tx2gene = tx2gene)

        ########################################################################

        read.counts=as.data.frame(txi.kallisto$counts)

        tpm=as.data.frame(txi.kallisto$abundance)

        efflen=as.data.frame(txi.kallisto$length)

        ########################################################################

        ## normalization

        norm.data=normalization(tpm)
        tpm.norm=norm.data[["expdata.norm"]]
        rownames(tpm.norm)=rownames(tpm)

        hk.genes=norm.data[["hk.genes"]]

        ########################################################################

        ## write output

        writeLines(hk.genes, con=paste(pathExpression, annot,  "/HousekeepingGenes.txt",sep=""))

        write.table(efflen, file=paste(pathExpression, annot,  "/AllSamples_KallistoEffectiveLength.txt",sep=""), row.names=T, col.names=T, quote=F, sep="\t")

        write.table(read.counts, file=paste(pathExpression, annot, "/AllSamples_KallistoEstimatedCounts.txt",sep=""), row.names=T, col.names=T, quote=F, sep="\t")

        write.table(tpm, file=paste(pathExpression, annot, "/AllSamples_KallistoRawTPM.txt",sep=""), row.names=T, col.names=T, quote=F, sep="\t")

        write.table(tpm.norm, file=paste(pathExpression, annot,  "/AllSamples_KallistoNormalizedTPM.txt",sep=""), row.names=T, col.names=T, quote=F, sep="\t")

        ########################################################################
    }
}
########################################################################
