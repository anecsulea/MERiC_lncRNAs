#!/bin/bash

######################################################################

pathLiftOver=../../data/liftOver_files
pathResults=../../results/CNA_analysis

######################################################################

for dataset in CRE AllExonsv6
do
    liftOver ${pathResults}/${dataset}_ampdel_hg19.bed ${pathLiftOver}/hg19ToHg38.over.chain.gz  ${pathResults}/${dataset}_ampdel_hg38.bed ${pathResults}/${dataset}_ampdel_hg38.unmapped
done

######################################################################
