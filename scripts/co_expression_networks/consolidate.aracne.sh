#!/bin/bash

export exp_dataset=$1
export TF_dataset=$2
export path=$3
export pathTools=$4

####################################################################################

export pathResults=${path}/results/co_expression_networks/${exp_dataset}/${TF_dataset}

####################################################################################

java -Xmx64G -jar ${pathTools}/ARACNe-AP/dist/aracne.jar -o ${pathResults} --consolidate

####################################################################################
