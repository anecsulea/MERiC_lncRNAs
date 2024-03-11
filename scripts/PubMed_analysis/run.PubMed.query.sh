#!/bin/bash

export path=$1

###########################################################################

export pathResults=${path}/results/PubMed_search

###########################################################################

export today=`date +"%d_%m_%Y"`

esearch -db pubmed -query "hepatocellular carcinoma"[Title] | efetch -format medline > ${pathResults}/PubMed_hepatocellular_carcinoma_Title_${today}.medline 

###########################################################################

## there should be 76,491 PubMed entries
