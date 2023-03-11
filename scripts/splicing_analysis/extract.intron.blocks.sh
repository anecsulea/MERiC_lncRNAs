#!/bin/bash

export annot=$1
export path=$2

#####################################################################

export pathEnsembl=${path}/data/ensembl_annotations
export pathGeneOverlaps=${path}/results/gene_overlaps
export pathResults=${path}/results/splicing_analysis
export pathScripts=${path}/scripts/splicing_analysis

export release=109

#####################################################################

if [ -e ${pathResults} ]; then
    echo "path output exists"
else
    mkdir -p ${pathResults}/
fi

#####################################################################

perl ${pathScripts}/extract.intron.blocks.pl --pathExonBlocks=${pathGeneOverlaps}/ExonBlocks_${annot}.txt --pathGeneInfo=${pathEnsembl}/GeneInfo_Ensembl${release}.txt --pathSynonyms=NA --biotypes=all --pathOutputIntronBlocks=${pathResults}/IntronBlocks_${annot}.txt 

#####################################################################
