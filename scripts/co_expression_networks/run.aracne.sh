#!/bin/bash

export exp_dataset=$1
export TF_dataset=$2
export path=$3
export pathTools=$4
export threads=$5

####################################################################################

export pathResults=${path}/results/co_expression_networks/${exp_dataset}

####################################################################################

for i in {1..100}
do
    echo $i

    java -Xmx5G -jar ${pathTools}/ARACNe-AP/dist/aracne.jar -e ${pathResults}/expression_data.txt -o ${pathResults}/${TF_dataset} --tfs ${pathResults}/${TF_dataset}.txt --pvalue 1E-8 --seed ${i} --threads=${threads}

done

####################################################################################
