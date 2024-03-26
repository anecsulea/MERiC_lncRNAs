#!/bin/bash

export annot=$1
export phast=$2
export path=$3
export scheduler=$4

#####################################################################

export pathEnsembl=${path}/data/ensembl_annotations
export pathGeneOverlaps=${path}/results/gene_overlaps
export pathResults=${path}/results/sequence_conservation/${phast}
export pathPhastCons=${path}/data/PhastCons/${phast}
export pathScripts=${path}/scripts/sequence_conservation

export release=109

#####################################################################

if [ ${annot} = "EnsemblNoOverlaps" ]; then
    export pathExons=${pathGeneOverlaps}
    export suffixExons=ExonBlocks_ExcludingOverlapOtherGenes_AllTranscripts_Ensembl${release}
fi

#####################################################################

if [ -e ${pathResults} ]; then
    echo "dir output already there"
else
    mkdir ${pathResults}
fi

#####################################################################

export suffix=.phastCons${phast}.wigFix.gz

#####################################################################

for chr in {1..22} X Y
do
    export pathPhast=${pathPhastCons}/chr${chr}${suffix}

    if [ -e ${pathPhast} ]; then

	if [ ${scheduler} = "slurm" ]; then
	    echo "#!/bin/bash " > ${pathScripts}/bsub_script
	    echo "#SBATCH --job-name=phast_${chr}" >>  ${pathScripts}/bsub_script
	    echo "#SBATCH --partition=normal" >>  ${pathScripts}/bsub_script
	    echo "#SBATCH --output=${pathScripts}/std_out_phast_${chr}" >>  ${pathScripts}/bsub_script
	    echo "#SBATCH --error=${pathScripts}/std_err_phast_${chr}" >>  ${pathScripts}/bsub_script
	    echo "#SBATCH --cpus-per-task=1" >>  ${pathScripts}/bsub_script ## ${ncores} CPU
	    echo "#SBATCH --time=24:00:00" >>  ${pathScripts}/bsub_script ## 24 hours

	    if [ ${chr} = Y ]||[ ${chr} = 19 ]||[ ${chr} = 20 ]; then
		echo "#SBATCH --mem=10G" >>  ${pathScripts}/bsub_script ## 10g per CPU
	    else
		echo "#SBATCH --mem=30G" >>  ${pathScripts}/bsub_script ## 30g per CPU
	    fi

	    echo "perl ${pathScripts}/compute.phastcons.exons.pl --pathExonBlocks=${pathExons}/${suffixExons}.txt --pathPhastCons=${pathPhast} --chr=${chr} --pathOutputExons=${pathResults}/PhastCons_ExonAverage_chr${chr}_${annot}.txt --pathOutputGenes=${pathResults}/PhastCons_GeneAverage_chr${chr}_${annot}.txt " >> bsub_script
	    
	    sbatch ${pathScripts}/bsub_script
	    
	fi
	
	if [ ${scheduler} = "none" ]; then
	    perl ${pathScripts}/compute.phastcons.exons.pl --pathExonBlocks=${pathExons}/${suffixExons}.txt --pathPhastCons=${pathPhast} --chr=${chr} --pathOutputExons=${pathResults}/PhastCons_ExonAverage_chr${chr}_${annot}.txt --pathOutputGenes=${pathResults}/PhastCons_GeneAverage_chr${chr}_${annot}.txt
	fi
    fi
done

#####################################################################
