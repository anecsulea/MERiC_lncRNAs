#!/bin/bash

export annot=$1
export path=$2

#########################################################################

export pathEnsembl=${path}/data/ensembl_annotations
export pathResults=${path}/results/gene_overlaps
export pathScripts=${path}/scripts/gene_overlaps

export release=109

#########################################################################

if [ ${annot} = "Ensembl" ]; then
    export pathGTF=${pathEnsembl}/Allranscripts_Ensembl${release}.gtf
    export suffix=AllTranscripts_Ensembl${release}
fi

#########################################################################

perl ${pathScripts}/distance.tss.pl --pathGTF=${pathGTF} --maxDistance=1000 --pathOutput=${pathResults}/DistanceAntisenseTSS_MaxDist1kb_${suffix}.txt

perl ${pathScripts}/distance.tss.pl --pathGTF=${pathGTF} --maxDistance=5000 --pathOutput=${pathResults}/DistanceAntisenseTSS_MaxDist5kb_${suffix}.txt

#########################################################################
