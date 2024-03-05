#!/bin/bash

export path=$1

#####################################################################

export pathResults=${path}/results/PubMed_search
export pathEnsembl=${path}/data/ensembl_annotations
export pathScripts=${path}/scripts/PubMed_analysis

export release=109

#####################################################################

export pathsPubMed=""

for file in `ls ${pathResults} | grep ^hepatocellular_carcinoma_Title`
do
    export pathsPubMed=${pathResults}/${file},${pathsPubMed}
done

#####################################################################

perl ${pathScripts}/parse.PubMed.results.pl --pathPubMedResults=${pathsPubMed} --pathGeneInfo=${pathEnsembl}/GeneInfo_Ensembl${release}.txt --pathGeneNames=${pathEnsembl}/GeneInfo_Ensembl${release}.txt --pathSynonyms=${pathEnsembl}/GeneNameSynonyms_Ensembl${release}.txt --pathForbiddenGenes=${pathResults}/forbidden_genes.txt  --pathOutput=${pathResults}/formatted_results_hepatocellular_carcinoma_Title.txt

#####################################################################
