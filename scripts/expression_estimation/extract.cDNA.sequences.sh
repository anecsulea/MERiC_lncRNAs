#!/bin/bash

export annot=$1

###########################################################################

export pathEnsembl=../../data/ensembl_annotations
export pathGenome=../../data/genome_sequence
export pathResults=../../results/kallisto_indexes
export pathScripts=../../scripts/expression_estimation

export release=109
export version=hg38

##############################################################

if [ -e ${pathResults} ]; then
    echo "output dir exists"
else
    mkdir -p ${pathResults}
fi

##############################################################

if [ ${annot} = "AllTranscripts" ]; then
    export pathGTF=${pathEnsembl}
    export suffix=AllTranscripts_Ensembl${release}
fi

##############################################################

export haplo=""

for chr in `cut -f 3 ${pathEnsembl}/GeneInfo_Ensembl${release}.txt | sort -u | grep CHR`
do
    export haplo=${chr},${haplo}
done

##############################################################

perl ${pathScripts}/extract.cDNA.sequences.pl --pathAnnotGTF=${pathGTF}/${suffix}.gtf.gz --forbiddenChromo=MT,${haplo} --pathGeneInfo=${pathEnsembl}/GeneInfo_Ensembl${release}.txt --forbiddenBiotypes=Mt_rRNA,rRNA,rRNA_pseudogene --pathGenomeSequence=${pathGenome}/genome_ensembl${release}_${version}.fa.gz --pathOutput=${pathResults}/${suffix}_noMT_norRNA_nohaplo.fa

###############################################################
