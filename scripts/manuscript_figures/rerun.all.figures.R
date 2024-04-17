############################################################

objects=ls()

files=c(paste("figure", 1:5, ".R",sep=""), paste("suppfigure", 1:8, ".R",sep=""))

############################################################

for(file in files){

    clear=setdiff(ls(), c("objects", "file", "files"))

    rm(list=clear)

    print(file)

    source(file)
}

############################################################
