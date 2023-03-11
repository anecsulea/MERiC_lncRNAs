#!/bin/bash

export path=$1

########################################################################

export pathMiTranscriptome=${path}/data/MiTranscriptome/hg38
export pathCHESS=${path}/data/CHESS
export pathEnsembl=${path}/data/ensembl_annotations
export pathHisat=${path}/data/hisat2_index

########################################################################

hisat2-build --seed 19 -p 32 --ss ${pathHisat}/splicesites.txt --exon ${pathHisat}/exons.txt ${pathHisat}/genome_sequence.fa ${pathHisat}/genome_sequence

########################################################################
