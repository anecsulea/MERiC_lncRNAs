############################################################

objects=ls()

files=c("data.gene.info.R", "data.PubMed.analysis.R",  "data.diffexp.R",  "data.expression.levels.R", "data.sequence.conservation.R", "data.gene.overlaps.R")

############################################################

for(file in files){

    clear=setdiff(ls(), c("objects", "file", "files"))

    rm(list=clear)

    print(file)

    source(file)
}

############################################################
