#!/bin/bash

export annot=$1
export windowsize=$2

###########################################################################

export pathAnnot=../../data/ensembl_annotations
export pathPCHiC=../../data/PCHiC_interactions
export pathResults=../../results/PCHiC_interactions
export pathScripts=../../scripts/PCHiC_analysis

###########################################################################

perl ${pathScripts}/extract.promoter.promoter.interactions.pl --pathAnnotGTF=${pathAnnot}/${annot}.gtf.gz  --pathGeneInfo=${pathEnsembl}/GeneInfo_Ensembl${release}.txt --pathInteractions=${pathPCHiC}/all_interactions.txt --windowSize=${windowsize} --interactionType=unbaited --pathOutput=${pathResults}/PromoterPromoterInteractions_${annot}_WindowSize${windowsize}_unbaited.txt

###########################################################################
