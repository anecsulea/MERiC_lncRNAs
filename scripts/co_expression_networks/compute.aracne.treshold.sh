#!/bin/bash

export exp_dataset=$1
export TF_dataset=$2
export path=$3
export pathTools=$4

####################################################################################

export pathResults=${path}/results/co_expression_networks/${exp_dataset}

####################################################################################

if [ -e ${pathResults}/${TF_dataset} ]; then
    echo "outdir already there"
else
    mkdir -p ${pathResults}/${TF_dataset}
fi

####################################################################################

java -Xmx5G -jar ${pathTools}/ARACNe-AP/dist/aracne.jar -e ${pathResults}/expression_data.txt  -o ${pathResults}/${TF_dataset} --tfs ${pathResults}/${TF_dataset}.txt --pvalue 1E-8 --seed 1 --calculateThreshold

####################################################################################
