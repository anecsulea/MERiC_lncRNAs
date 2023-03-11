######################################################################

pathData="../../data/CNA/2022/"
pathResults="../../results/CNA_analysis/"

######################################################################

for(dataset in c("CRE", "AllExonsv6")){
    orig=read.table(paste(pathData, "all_",dataset,".geneCN.ampdel.txt",sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
    orig=orig[which(!is.na(orig$seg_start)),]
    orig$id=paste(orig$seg_chr,":", orig$seg_start, "-", orig$seg_end,sep="")
    orig$sample=unlist(lapply(orig$seg_sample, function(x) unlist(strsplit(x, split="_GL"))[1]))
    orig$method=unlist(lapply(orig$seg_sample, function(x) unlist(strsplit(x, split="_GL_"))[2]))
    orig$tumor_sample=unlist(lapply(orig$sample, function(x) unlist(strsplit(x, split="_"))[1]))
    orig$normal_sample=unlist(lapply(orig$sample, function(x) unlist(strsplit(x, split="_"))[2]))

    print(paste(nrow(orig), "lines originally"))

    lifted=read.table(paste(pathResults,dataset, "_ampdel_hg38.bed",sep=""),h=F, stringsAsFactors=F)
    colnames(lifted)=c("chr", "start", "end", "id", "strand")
    lifted$chr=unlist(lapply(lifted$chr, function(x) substr(x, 4, nchar(x))))

    lifted$newcoords=paste(lifted$chr, ":", lifted$start, "-", lifted$end,sep="")
    lifted$completeid=paste(lifted$id, lifted$newcoords)
    lifted=lifted[which(!duplicated(lifted$completeid)),]

    ## keep only largest lifted segment
    lifted$size=lifted$end-lifted$start
    lifted=lifted[order(lifted$size, decreasing=T),]
    lifted=lifted[order(lifted$id),]
    lifted=lifted[which(!duplicated(lifted$id)),]

    ## old chr = new chr
    lifted$oldchr=unlist(lapply(lifted$id, function(x) unlist(strsplit(x, split=":"))[1]))
    lifted=lifted[which(lifted$oldchr==lifted$chr),]

    ## filted on size ratio
    lifted$oldcoords=unlist(lapply(lifted$id, function(x) unlist(strsplit(x, split=":"))[2]))
    lifted$oldstart=as.numeric(unlist(lapply(lifted$oldcoords, function(x) unlist(strsplit(x, split="-"))[1])))
    lifted$oldend=as.numeric(unlist(lapply(lifted$oldcoords, function(x) unlist(strsplit(x, split="-"))[2])))
    lifted$oldsize=lifted$oldend-lifted$oldstart
    lifted$sizeratio=lifted$oldsize/lifted$size
    lifted=lifted[which(lifted$sizeratio>=0.9 & lifted$sizeratio<=1.1),]

    ##  combine results add sample etc. info
    rownames(lifted)=lifted$id

    orig=orig[which(orig$id%in%lifted$id),]
    results=orig[,c("seg_chr", "seg_start", "seg_end", "tumor_sample", "normal_sample", "seg_type", "method")]
    results$chr_hg38=lifted[orig$id, "chr"]
    results$start_hg38=lifted[orig$id, "start"]
    results$end_hg38=lifted[orig$id, "end"]

    print(paste(nrow(results), "lines kept"))

    write.table(results, file=paste(pathResults, "all_", dataset,".geneCN.ampdel.hg38.txt",sep=""), sep="\t", row.names=F, col.names=T, quote=F)

}

######################################################################
