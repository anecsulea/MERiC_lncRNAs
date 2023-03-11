#!/bin/bash

export path=$1

########################################################################

export pathMiTranscriptome=${path}/data/MiTranscriptome/hg38
export pathCHESS=${path}/data/CHESS
export pathEnsembl=${path}/data/ensembl_annotations

########################################################################

hisat2_extract_exons.py ${pathMiTranscriptome}/mitranscriptome.v2.exons.filtered.gtf > ${pathMiTranscriptome}/mitranscriptome.v2.exons.hisat2.txt

hisat2_extract_splice_sites.py ${pathMiTranscriptome}/mitranscriptome.v2.exons.filtered.gtf > ${pathMiTranscriptome}/mitranscriptome.v2.splicesites.hisat2.txt

########################################################################

hisat2_extract_exons.py ${pathCHESS}/chess2.2.filtered.gtf > ${pathCHESS}/chess2.2.exons.hisat2.txt

hisat2_extract_splice_sites.py ${pathCHESS}/chess2.2.filtered.gtf > ${pathCHESS}/chess2.2.splicesites.hisat2.txt

########################################################################

hisat2_extract_exons.py ${pathEnsembl}/AllTranscripts_Ensembl97.gtf > ${pathEnsembl}/AllTranscripts_Ensembl97.exons.hisat2.txt

hisat2_extract_splice_sites.py ${pathEnsembl}/AllTranscripts_Ensembl97.gtf > ${pathEnsembl}/AllTranscripts_Ensembl97.splicesites.hisat2.txt

########################################################################
