#############################################################################

pathExpression="../../results/expression_estimation/"
pathDifferentialExpression="../../results/differential_expression/"
pathAnnot="../../data/ensembl_annotations/"
pathSampleInfo="../../results/sample_info/"
pathPubMed="../../results/PubMed_search/"
pathResults="../../results/figures/main_figures/"

annot="AllTranscripts_Ensembl109_noMT_norRNA_nohaplo"
release=109

#############################################################################

if(prepare){
  
  ###################################################################################
  ## nb cited
  
  nb.cit=read.table(paste(pathPubMed, "number_of_citations_hepatocellular_carcinoma_Title.txt",sep=""),h=T, stringsAsFactors=F)

  top.cited.pc=nb.cit$GeneID[which(nb.cit$GeneType=="protein_coding")][1:10]
  highly.cited.lnc=nb.cit$GeneID[which(nb.cit$GeneType=="lncRNA")][1:10]

  ###################################################################################
  
  sampleinfo=read.table(paste(pathSampleInfo, "AnalyzedSamples.txt", sep=""), h=T, stringsAsFactors=F,sep="\t")
  
  sampleinfo$NAFLD=rep("no", nrow(sampleinfo))
  sampleinfo$NAFLD[grep("ALD", sampleinfo$Disease)]="yes"
  
  sampleinfo$ALD=rep("no", nrow(sampleinfo))
  sampleinfo$ALD[grep("ALD", sampleinfo$Disease)]="yes"
  
  sampleinfo$HepB=rep("no", nrow(sampleinfo))
  sampleinfo$HepB[grep("Hepatitis B", sampleinfo$Disease)]="yes"
  
  sampleinfo$HepC=rep("no", nrow(sampleinfo))
  sampleinfo$HepC[grep("Hepatitis C", sampleinfo$Disease)]="yes"
  
  ## sample info, only for tumor samples

  tumor.samples=sampleinfo[which(sampleinfo$TissueType=="Tumor"),]
  liver.samples=sampleinfo[which(sampleinfo$TissueType=="Liver"),]
  
  rownames(tumor.samples)=as.character(tumor.samples$PatientID)
    
  ## we order them by all the factors
  tumor.samples=tumor.samples[order(tumor.samples$Sex),]
  tumor.samples=tumor.samples[order(tumor.samples$NAFLD),]
  tumor.samples=tumor.samples[order(tumor.samples$ALD),]
  tumor.samples=tumor.samples[order(tumor.samples$HepB),]
  tumor.samples=tumor.samples[order(tumor.samples$HepC),]
  tumor.samples=tumor.samples[order(tumor.samples$Cirrhosis),]
  tumor.samples=tumor.samples[order(tumor.samples$EdmondsonGrade),]

 
########################################################################
  
  de.tl=read.table(paste(pathDifferentialExpression,annot, "/DifferentialExpression_Tumor_vs_Liver.txt", sep=""), h=T, stringsAsFactors=F, sep="\t")
  
  de.grades=read.table(paste(pathDifferentialExpression,annot, "/DifferentialExpression_EdmondsonGrade_34_vs_12.txt", sep=""), h=T, stringsAsFactors=F, sep="\t")

  diffexp=list("Tumor_vs_NonTumor"=de.tl, "EdmondsonGrade"=de.grades)

  maxFDR=0.05
  minFC=1.5
  
  ## DE table

  de.table=matrix(rep(NA, length(highly.cited.lnc)*2), nrow=length(highly.cited.lnc))
  
  rownames(de.table)=highly.cited.lnc
  colnames(de.table)=c("Tumor_vs_NonTumor", "EdmondsonGrade")

  for(test in colnames(de.table)){
    this.de=diffexp[[test]]
    tested=rownames(this.de)
    signif.up=rownames(this.de)[which(this.de$padj<maxFDR & this.de$log2FoldChange>log2(minFC))]
    signif.down=rownames(this.de)[which(this.de$padj<maxFDR & this.de$log2FoldChange<log2(1/minFC) )]
    
    relax.up=rownames(this.de)[which(this.de$padj<maxFDR & this.de$log2FoldChange>0)]
    relax.down=rownames(this.de)[which(this.de$padj<maxFDR & this.de$log2FoldChange<0 )]

    relax.up=setdiff(relax.up, signif.up)
    relax.down=setdiff(relax.down, signif.down)

    if(any(highly.cited.lnc%in%tested)){
      de.table[intersect(highly.cited.lnc, tested), test]=0
    }
    
    if(any(highly.cited.lnc%in%signif.up)){
      de.table[intersect(highly.cited.lnc, signif.up), test]=1
    } 
    if(any(highly.cited.lnc%in%relax.up)){
      de.table[intersect(highly.cited.lnc, relax.up), test]=0.5
    }
    
    if(any(highly.cited.lnc%in%signif.down)){
      de.table[intersect(highly.cited.lnc, signif.down), test]=-1
    }
    
    if(any(highly.cited.lnc%in%relax.down)){
      de.table[intersect(highly.cited.lnc, relax.down), test]=-0.5
    }
  }

  ## gene exp
  tpm=read.table(paste(pathExpression, annot, "/AllSamples_KallistoRawTPM.txt",sep=""), h=T, stringsAsFactors=F)
  log2tpm=log2(tpm+1)

  
  col.Edmondson=c("gold", "orange", "darkorange3", "red")
  prepare=FALSE
}

#############################################################################

## 1 column width 87 mm = 3.34 in
## 2 columns width 180 mm = 7.04 in
## max height: 9 in, including legend

#############################################################################

pdf(file=paste(pathResults, "Figure3.pdf", sep=""), width=6.85, height=6.95)

m=matrix(rep(NA, 32*19), nrow=32)

## first line: spacer
for(i in 1){
  m[i,]=rep(1, 19)
}

## individual lines:  Edmondson grade, cirrhosis, hepC, hepB, ALD, NAFLD, sex 

for(i in 2:8){
  m[i,]=c(rep(2*(i-1), 6),  rep(2*(i-1)+1, 13))
}

## spacer plot

for(i in 9:10){
  m[i,]=rep(16, 19)
}

## expression of highly cited lncRNAs


## a*11+b=17
## a*13+b=20

a=(20-17)/2
b=17-11*a

for(i in seq(from=11, by=2, length.out=length(highly.cited.lnc))){
  j=a*i+b
  m[i,]=c(rep(j,5), rep(j+1,1),  rep(j+2,13))
  m[i+1,]=c(rep(j,5), rep(j+1,1),  rep(j+2,13))
}

## legend 

for(i in 31:32){
  M=max(m, na.rm=T)
  m[i,]=rep(M+1, 19)
}

layout(m)

#############################################################################

## legend

par(mar=c(0,0,0,0))
plot(1, type="n", xlab="", ylab="", axes=F, xlim=c(0,1))

mtext("A", side=3, at=-0.025, line=-1.2, cex=0.95, font=2)

#############################################################################

labels=c("Edmondson grade", "Cirrhosis", "Hepatitis C", "Hepatitis B", "ALD", "NAFLD", "Sex")
names(labels)=c("EdmondsonGrade",  "Cirrhosis", "HepC", "HepB", "ALD", "NAFLD","Sex")

for(column in c("EdmondsonGrade", "Cirrhosis", "HepC", "HepB", "ALD", "NAFLD", "Sex")){

  if(column=="BiopsyAge"){
    fac.age=cut(tumor.samples[,column], breaks=c(0, 50, 60, 70, 80, 90), include.lowest=T)
    col.age=gray.colors(n=length(levels(fac.age)), rev=TRUE)
    names(col.age)=levels(fac.age)
    
    this.colors=col.age[as.character(fac.age)]
    
    levels=c("<50", "50-60", "60-70", "70-80", ">80")
  } else{
    if(column=="EdmondsonGrade"){
      col.factor=col.Edmondson
      this.colors=col.Edmondson[tumor.samples[,"EdmondsonGrade"]]
      levels=1:4
    } else{
    
      col.factor=gray.colors(n=length(levels(as.factor(tumor.samples[,column]))), rev=TRUE)
      names(col.factor)=levels(as.factor(tumor.samples[,column]))
      
      this.colors=col.factor[as.character(tumor.samples[,column])]
      
      levels=names(col.factor)
    }
  }
  
  label=labels[column]
   
  ## legend for this plot
  
  par(mar=c(0,0,0,0))
  plot(1, type="n", xlab="", ylab="", axes=F, ylim=c(0,1), xlim=c(0, 1), xaxs="i", yaxs="i")

  ## mtext(label, side=3, cex=0.5, line=0, at=0)

  if(column=="BiopsyAge"){
    text(label, x=0.1, y=0.5, adj=c(0, 0.5), cex=0.95, xpd=NA)
    legend("topleft", legend=levels[1], fill=col.age[1], border=col.age[1], bty="n", horiz=T, xpd=NA, cex=0.85, inset=c(0.2, -0.3))
    legend("topleft", legend=levels[2:4], fill=col.age[2:4], border=col.age[2:4], bty="n", horiz=T, xpd=NA, cex=0.85, inset=c(0.4, -0.3))
    legend("topleft", legend=levels[5], fill=col.age[5], border=col.age[5], bty="n", horiz=T, xpd=NA, cex=0.85, inset=c(0.2, 0.55))
  } else{
    if(column=="EdmondsonGrade"){
      legend("topleft", legend=levels, fill=col.factor, border=col.factor, bty="n", horiz=T, xpd=NA, cex=0.85, inset=c(0.55, -0.32))
      text(label, x=0.1, y=0.55, adj=c(0, 0.5), cex=0.95, xpd=NA)
    } else{
      legend("topleft", legend=levels[2:1], fill=col.factor[2:1], border=col.factor[2:1], bty="n", horiz=T, xpd=NA, cex=0.85, inset=c(0.55, -0.32))
      text(label, x=0.1, y=0.55, adj=c(0, 0.5), cex=0.95, xpd=NA)
     }
  }

  ## actual plot
  
  par(mar=c(0.1, 2.8, 0.15, 3))
  image(cbind(1:dim(tumor.samples)[1]), col=this.colors, axes=F, xaxs="i", yaxs="i")
  
}

#############################################################################

## spacer plot

plot(1, type="n", xlab="", ylab="", axes=F, xlim=c(0,1))
mtext("B", side=3, at=-0.087, line=-1.25, cex=0.95, font=2)
mtext("C", side=3, at=0.25, line=-1.25, cex=0.95, font=2)
mtext("D", side=3, at=0.31, line=-1.25, cex=0.95, font=2)
mtext("E", side=3, at=1.055, line=-1.25, cex=0.95, font=2)

#############################################################################

## selected lncRNAs

for(id in highly.cited.lnc){
  ## gene name
  name=geneinfo[id,"Gene.name"]

  print(name)

  ## gene expression values for sample categories
  this.exp.tumor=as.numeric(log2tpm[id, tumor.samples$BiopsyID])
  this.exp.liver=as.numeric(log2tpm[id, liver.samples$BiopsyID])
    
  ## difference between tumor and liver

  par(mar=c(1.5, 4, 0.1, 0.15))
  boxplot(this.exp.liver,this.exp.tumor,  col="white",  border=c("gray40", "red"), axes=F, outline=F, boxwex=0.5, horizontal=T)
  
  par(tck=-0.15)  
  axis(side=1, mgp=c(3, 0.1, 0), cex.axis=0.7)
  box()
  mtext(name, side=2, las=2, cex=0.6, line=0.15, font=3, adj=1) 
  
  ## differential expression

  par(mar=c(0.5,1.65,0.1,0))
  xlim=c(0,1)
  ylim=c(0,1)
  plot(1, type="n", xlab="", ylab="",axes=F, xlim=xlim, ylim=ylim, xaxs="i", yaxs="i")

  ## tumor vs normal liver

  x2=0.5
 
  ## tumor vs non-tumor, paired samples
  
  if(de.table[id,"Tumor_vs_NonTumor"]==1){
    arrows(x2, 0.1, x2, 0.9, length=0.035, lwd=1.15, col="black")
  }
  if(de.table[id,"Tumor_vs_NonTumor"]==0.5){
    arrows(x2, 0.1, x2, 0.9, length=0.035, lwd=1.15, col="gray50")
  }

  if(de.table[id,"Tumor_vs_NonTumor"]==-1){
    arrows(x2, 0.9, x2, 0.1, length=0.035, lwd=1.15, col="black")
  }
   if(de.table[id,"Tumor_vs_NonTumor"]==-0.5){
    arrows(x2, 0.9, x2, 0.1, length=0.035, lwd=1.15, col="gray50")
  }
    
  ## plot for tumors only
  
  ylim=c(0, max(c(log2(this.exp.tumor+1),1)))
  
  xlim=c(0.25, length(this.exp.tumor)+0.5)
  xwidth=0.25
  
  par(mar=c(0.5, 2.8, 0.1, 3))
  plot(1, type="n", xlab="", ylab="",axes=F, xlim=xlim, ylim=ylim, xaxs="i")

  xpos=1:length(this.exp.tumor)
  this.values=log2(this.exp.tumor+1)
  this.colors=col.Edmondson[tumor.samples$EdmondsonGrade]
  rect(xpos-xwidth, 0, xpos+xwidth, this.values, col=this.colors, border=this.colors)

  ## abline(v=max(xpos), xpd=NA, lwd=0.2) ## to check alignment
  ## abline(v=min(xpos), xpd=NA) ## to check alignment

  maxval=ylim[2]
  par(tck=-0.25)
  axis(side=1, at=range(xlim)+c(-1,1), labels=rep("",2))
  axis(side=2, at=c(0, ylim[2]), labels=rep("",2), cex.axis=0.5, mgp=c(3, 0.65, 0))
  mtext(c(0, round(ylim[2], digits=1)), at=c(0, ylim[2]), las=2, cex=0.5, line=0.65, adj=1, side=2)

  ##text(name, x=xlim[2], adj=1, y=ylim[2]+diff(ylim)/10, cex=0.85, xpd=NA)
  
  ## DE, Edmondson grade
  
  x3=xlim[2]+diff(xlim)/20
  y1=diff(ylim)/10
  y2=ylim[2]-diff(ylim)/10
   
  if(de.table[id,"EdmondsonGrade"]==1){
    arrows(x3, y1, x3, y2, length=0.035, lwd=1.15, xpd=NA, col="black")
  }
  if(de.table[id,"EdmondsonGrade"]==0.5){
    arrows(x3, y1, x3, y2, length=0.035, lwd=1.15, xpd=NA, col="gray50")
  }

  if(de.table[id,"EdmondsonGrade"]==-1){
    arrows(x3, y2, x3, y1, length=0.035, lwd=1.15, xpd=NA, col="black")
  }
  if(de.table[id,"EdmondsonGrade"]==-0.5){
    arrows(x3, y2, x3, y1, length=0.035, lwd=1.15, xpd=NA, col="gray50")
  }
}

#############################################################################

par(mar=c(0, 0, 0, 0))
plot(1, type="n", xlab="", ylab="", axes=F, xlim=c(0,1), ylim=c(0,1))

text("1 & 2", x=1.01, y=1.1, cex=0.9, xpd=NA)
text("vs.", x=1.01, y=0.38, cex=0.9, xpd=NA, font=3)
text("3 & 4", x=1.01, y=-0.36, cex=0.9, xpd=NA)

text("log2(TPM+1)", x=0.15, y=0.3, cex=0.9, xpd=NA)

legend("topleft", lty=1, col=c("gray40", "red"), legend=c("tumor", "liver"), cex=0.85, bty="n", xpd=NA)

## legend for colors

text("tumor", x=0.288, y=0.45, cex=0.9, xpd=NA)
text("vs.", x=0.288, y=-0.2, cex=0.9, xpd=NA, font=3)
text("liver", x=0.288, y=-0.8, cex=0.9, xpd=NA)

legend("topleft", legend=as.character(1:4), fill=col.Edmondson, border=col.Edmondson, bty="n", cex=0.85, inset=c(0.77, 0), horiz=T, xpd=NA)


text("log2(TPM+1)", x=0.37, y=0.3, cex=0.9, xpd=NA)
text("Edmondson grade", x=0.72, y=0.3, cex=0.9, xpd=NA)

## segments(0, 0.1, 0, 0.7, col="red", xpd=NA)
## text("median", x=0, y=-0.59, cex=0.9, xpd=NA)

#############################################################################

dev.off()

#############################################################################
