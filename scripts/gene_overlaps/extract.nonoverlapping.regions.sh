#!/bin/bash

export annot=$1
export path=$2

##############################################################

export pathEnsembl=${path}/data/ensembl_annotations
export pathResults=${path}/results/gene_overlaps
export pathScripts=${path}/scripts/gene_overlaps

export release=109

##############################################################

if [ ${annot} = "Ensembl" ]; then
    export pathAnnot=${pathEnsembl}
    export prefixAnnot=AllTranscripts_Ensembl${release}
fi

##############################################################

perl ${pathScripts}/extract.nonoverlapping.regions.pl --pathAnnotGTF=${pathAnnot}/${prefixAnnot}.gtf --pathRepeats=NA --pathRetrogenes=NA --pathOutput=${pathResults}/ExonCoordinates_ExcludingOverlapOtherGenes_${prefixAnnot}.txt

##############################################################
