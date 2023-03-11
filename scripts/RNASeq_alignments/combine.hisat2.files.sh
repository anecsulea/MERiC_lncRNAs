#!/bin/bash

export path=$1

########################################################################

export pathMiTranscriptome=${path}/data/MiTranscriptome/hg38
export pathCHESS=${path}/data/CHESS
export pathEnsembl=${path}/data/ensembl_annotations
export pathHisat=${path}/data/hisat2_index

########################################################################

cat ${pathEnsembl}/AllTranscripts_Ensembl97.splicesites.hisat2.txt ${pathCHESS}/chess2.2.splicesites.hisat2.txt ${pathMiTranscriptome}/mitranscriptome.v2.splicesites.hisat2.txt | sort -u > ${pathHisat}/splicesites.txt

########################################################################

cat ${pathEnsembl}/AllTranscripts_Ensembl97.exons.hisat2.txt ${pathCHESS}/chess2.2.exons.hisat2.txt ${pathMiTranscriptome}/mitranscriptome.v2.exons.hisat2.txt | sort -u > ${pathHisat}/exons.txt

########################################################################
