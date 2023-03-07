#!/bin/bash

export annot=$1
export path=$2
export pathTools=$3

####################################################################################

export pathResults=${path}/results/kallisto_indexes

####################################################################################

singularity exec -B ${path} ${pathTools}/kallisto.sif kallisto index -i ${pathResults}/${annot} ${pathResults}/${annot}.fa

####################################################################################
