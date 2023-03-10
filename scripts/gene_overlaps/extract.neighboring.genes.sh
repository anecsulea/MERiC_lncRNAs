#!/bin/bash

export annot=$1
export path=$2

#########################################################################

export pathResults=${path}/results/gene_overlaps
export pathScripts=${path}/scripts/gene_overlaps

#########################################################################

for dist in 50 100 250 500
do
    perl ${pathScripts}/extract.neighboring.genes.pl --pathExonBlocks=${pathResults}/ExonBlocks_${annot}.txt --windowSize=${dist}000 --pathOutput=${pathResults}/NeighboringGenes_${dist}kb_${annot}.txt
done

#########################################################################

perl ${pathScripts}/extract.neighboring.genes.pl --pathExonBlocks=${pathResults}/ExonBlocks_${annot}.txt --windowSize=0 --pathOutput=${pathResults}/OverlappingGenes_${annot}.txt

#########################################################################







