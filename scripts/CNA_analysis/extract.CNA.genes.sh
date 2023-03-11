#!/bin/bash

export annot=$1
export path=$2

#########################################################################

export pathAnnot=${path}/results/gene_overlaps
export pathResults=${path}/results/CNA_analysis
export pathScripts=${path}/scripts/CNA_analysis

#########################################################################

perl ${pathScripts}/extract.CNA.genes.pl --pathExonBlocks=${pathAnnot}/ExonBlocks_${annot}.txt --pathsCNACoordinates=${pathResults}/all_AllExonv6.geneCN.ampdel.hg38.txt,${pathResults}/all_CRE.geneCN.ampdel.hg38.txt --pathOutput=${pathResults}/GeneOverlap_${annot}_CNA_AllExonv6_CRE.txt

#########################################################################
