#!/bin/bash

export annot=$1
export phast=$2
export path=$3

#####################################################################

export pathResults=${path}/results/sequence_conservation/${phast}
export pathScripts=${path}/scripts/sequence_conservation

#####################################################################

cp ${pathResults}/PhastCons_ExonAverage_chr1_${annot}.txt ${pathResults}/PhastCons_ExonAverage_${annot}.txt
cp ${pathResults}/PhastCons_GeneAverage_chr1_${annot}.txt ${pathResults}/PhastCons_GeneAverage_${annot}.txt

for chr in {2..22} X Y
do
    sed '1d' ${pathResults}/PhastCons_ExonAverage_chr${chr}_${annot}.txt >> ${pathResults}/PhastCons_ExonAverage_${annot}.txt
    sed '1d' ${pathResults}/PhastCons_GeneAverage_chr${chr}_${annot}.txt >> ${pathResults}/PhastCons_GeneAverage_${annot}.txt
done

#####################################################################
