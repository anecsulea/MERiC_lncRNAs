##################################################################

pathTables="../../data_for_publication/tables/"
pathEGA="../../data_for_publication/EGA_submission/"

##################################################################

samples=read.table(paste(pathTables, "SupplementaryTable5.tsv",sep=""), h=T, stringsAsFactors=F, sep="\t")

##################################################################

n=nrow(samples)

formatted=data.frame("alias"=paste("MERIC_adjacent_tissue", samples$biopsyID,sep="_"), title=samples$biopsyID, "description"=rep("RNA-seq from adjacent tissue biopsy",n), "biological_sex"=samples$sex, "subject_id"=samples$Patient_ID, "phenotype"=rep("hepatocellular carcinoma",n), "biosample_id"=rep(NA, n), "case_control"=rep("control", n), "organism_part"=rep("liver", n), "cell_line"=rep(NA, n))

formatted$biological_sex[which(formatted$biological_sex=="m")]="male"
formatted$biological_sex[which(formatted$biological_sex=="f")]="female"

write.table(formatted, file=paste(pathEGA, "EGA_samples.csv",sep=""), row.names=F, col.names=T, sep=",", quote=F)

##################################################################

