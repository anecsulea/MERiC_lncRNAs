#!/bin/bash

export annot=$1
export path=$2

#################################################################################

export pathGeneOverlaps=${path}/results/gene_overlaps
export pathScripts=${path}/scripts/gene_overlaps

################################################################################

perl make.exon.blocks.nooverlap.pl --pathExonCoords=${pathGeneOverlaps}/ExonCoordinates_ExcludingOverlapOtherGenes_${annot}.txt --collapseDistance=0 --pathOutput=${pathGeneOverlaps}/ExonBlocks_ExcludingOverlapOtherGenes_${annot}.txt

#################################################################################

    
