#!/bin/bash

export path=$1

#####################################################################

export pathAnnot=${path}/data/ensembl_annotations
export pahtResults=${path}/results/gene_overlaps
export pathScripts=${path}/scripts/gene_overlaps

export release=109

export suffix=AllTranscripts_Ensembl${release}

#####################################################################

if [ -e ${pathResults} ]; then
    echo "path output already there"
else
    mkdir -p ${pathResults}
fi

#####################################################################

perl ${pathScripts}/make.exon.blocks.gtf.pl --pathGTF=${pathAnnot}/${suffix}.gtf.gz --collapseDistance=0  --pathOutputExonBlocks=${pathResults}/ExonBlocks_${suffix}.txt

#################################################################################

    
