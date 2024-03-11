#!/bin/bash

export path=$1

#####################################################################

export pathResults=${path}/results/PubMed_search
export pathEnsembl=${path}/data/ensembl_annotations
export pathScripts=${path}/scripts/PubMed_analysis

export release=109

export date="11_03_2024"

#####################################################################

perl ${pathScripts}/parse.PubMed.results.medline.pl --pathPubMedResults=${pathResults}/PubMed_hepatocellular_carcinoma_Title_${date}.medline --pathGeneInfo=${pathEnsembl}/GeneInfo_Ensembl${release}.txt --pathGeneNames=${pathEnsembl}/GeneInfo_Ensembl${release}.txt --pathSynonyms=${pathEnsembl}/GeneNameSynonyms_Ensembl${release}.txt --pathForbiddenGenes=${pathResults}/forbidden_genes.txt  --pathOutput=${pathResults}/PubMed_hepatocellular_carcinoma_Title_${date}.txt

#####################################################################
