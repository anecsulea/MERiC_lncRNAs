########################################################################

pathResults="../../results/splicing_analysis/"
pathDocs="../../docs/"
pathFigures="../../results/figures/analyses/"

annot="AllTranscripts_Ensembl109"

########################################################################

tumor.samples=read.table(paste(pathDocs, "TumorSamples_Ng2022.txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
tumor.samples=tumor.samples[which(tumor.samples$tumor_biopsyID!=""),]

###########################################################################

liver.samples=read.table(paste(pathDocs, "LiverSamples_Ng2022.txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
liver.samples=liver.samples[which(liver.samples$biopsyID!="POOL"),]

###########################################################################

sample.info.t=data.frame("BiopsyID"=tumor.samples$tumor_biopsyID, "TissueType"=rep("Tumor", nrow(tumor.samples)), "Sex"=tumor.samples$sex, "AgeAtBiopsy"=tumor.samples$age_at_biopsy, "EdmondsonGrade"=tumor.samples$edmondson, "Cirrhosis"=tumor.samples$cirrhosis, "Diseases"=tumor.samples$underlying_liver_disease)

sample.info.l=data.frame("BiopsyID"=liver.samples$biopsyID, "TissueType"=rep("Liver", nrow(liver.samples)), "Sex"=liver.samples$sex, "AgeAtBiopsy"=liver.samples$age_at_biopsy, "EdmondsonGrade"=rep(NA, nrow(liver.samples)), "Cirrhosis"=rep(NA, nrow(liver.samples)), "Diseases"=rep(NA, nrow(liver.samples)))

sample.info=rbind(sample.info.t, sample.info.l)

samples=sample.info$BiopsyID

########################################################################

ir=read.table(paste(pathResults, "/AllSamples_IntronExonRatio_", annot, ".txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")

samples=intersect(colnames(ir), samples)

median.ir=apply(ir[,samples],2, median, na.rm=T)
names(median.ir)=samples

ts=intersect(tumor.samples$tumor_biopsyID, samples)
ls=intersect(liver.samples$biopsyID, samples)

########################################################################

density.tumor=density(median.ir[ts], bw=0.01)
density.liver=density(median.ir[ls], bw=0.01)

xlim=range(c(density.tumor$x, density.liver$x))
ylim=range(c(density.tumor$y, density.liver$y))

plot(density.tumor, col="red", xlim=xlim, ylim=ylim, main="", xlab="median intron/exon ratio", ylab="density")
lines(density.liver, col="black")

legend("topright", legend=c("liver", "tumor"), col=c("black", "red"), lty=1, inset=0.01)

########################################################################

pdf(file=paste(pathFigures, "IntronExonRatioHistogram.pdf",sep=""))

hist(median.ir[ts], col=alpha("red", 0.5), xlab="median intron/exon ratio", main="", breaks=seq(from=0, to=0.30, by=0.0125), xlim=c(0,0.25), ylab="number of samples")
hist(median.ir[ls], col=alpha("black", 0.5), add=T,breaks=seq(from=0, to=0.30, by=0.0125))

legend("topright", legend=c("liver", "tumor"), fill=c("black", "red"),inset=0.01)

dev.off()

########################################################################
