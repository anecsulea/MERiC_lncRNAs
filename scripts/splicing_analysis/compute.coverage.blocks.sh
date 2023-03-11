#!/bin/bash

export sample=$1
export annot=$2
export path=$3
export scheduler=$4

#####################################################################

export pathGeneOverlaps=${path}/results/gene_overlaps
export pathGenomeBrowser=${path}/results/genome_browser
export pathResults=${path}/results/splicing_analysis
export pathScripts=${path}/scripts/splicing_analysis

###########################################################################

if [ -e ${pathGenomeBrowser}/${sample}/coverage_unique_forward.bedGraph.gz ]&&[ -e ${pathGenomeBrowser}/${sample}/coverage_unique_reverse.bedGraph.gz ]; then
    
    if [ -e ${pathResults}/${sample} ]; then
	echo "path exists"
    else
	mkdir ${pathResults}/${sample}
    fi

    if [ -e ${pathResults}/${sample}/CoverageGenes_ExonBlocks_${annot}.txt ]; then
	echo "exons already done"
    else
	echo "#!/bin/bash " > ${pathScripts}/bsub_script

	if [ ${scheduler} = "slurm" ]; then
	    echo "#SBATCH --job-name=covex_${sample}" >>  ${pathScripts}/bsub_script
            echo "#SBATCH --output=${pathScripts}/std_out_covex_${sample}" >>  ${pathScripts}/bsub_script
	    echo "#SBATCH --error=${pathScripts}/std_err_covex_${sample}" >>  ${pathScripts}/bsub_script
	    echo "#SBATCH --cpus-per-task=1" >>  ${pathScripts}/bsub_script 
	    echo "#SBATCH --time=1:00:00" >>  ${pathScripts}/bsub_script 
	    echo "#SBATCH --mem=10G" >>  ${pathScripts}/bsub_script 
	fi
	
	echo "perl ${pathScripts}/compute.coverage.blocks.pl --pathExonBlocks=${pathGeneOverlaps}/ExonBlocks_${annot}.txt --pathCoverageForward=${pathGenomeBrowser}/${sample}/coverage_unique_forward.bedGraph.gz --pathCoverageReverse=${pathGenomeBrowser}/${sample}/coverage_unique_reverse.bedGraph.gz  --pathOutputGenes=${pathResults}/${sample}/CoverageGenes_ExonBlocks_${annot}.txt" >> ${pathScripts}/bsub_script

	if [ ${scheduler} = "slurm" ]; then
	    sbatch  ${pathScripts}/bsub_script
	else
	    chmod a+x ${pathScripts}/bsub_script
	    ${pathScripts}/bsub_script
	fi
    fi
    
    if [ -e ${pathResults}/${sample}/CoverageGenes_IntronBlocks_${annot}.txt ]; then
	echo "introns already done"
    else
	echo "#!/bin/bash " > ${pathScripts}/bsub_script

	if [ ${scheduler} = "slurm" ]; then
	    echo "#SBATCH --job-name=covint_${sample}" >>  ${pathScripts}/bsub_script
            echo "#SBATCH --output=${pathScripts}/std_out_covint_${sample}" >>  ${pathScripts}/bsub_script
	    echo "#SBATCH --error=${pathScripts}/std_err_covint_${sample}" >>  ${pathScripts}/bsub_script
	    echo "#SBATCH --cpus-per-task=1" >>  ${pathScripts}/bsub_script 
	    echo "#SBATCH --time=1:00:00" >>  ${pathScripts}/bsub_script 
	    echo "#SBATCH --mem=10G" >>  ${pathScripts}/bsub_script 
	fi
	
	echo "perl ${pathScripts}/compute.coverage.blocks.pl --pathExonBlocks=${pathResults}/IntronBlocks_${annot}.txt --pathCoverageForward=${pathGenomeBrowser}/${sample}/coverage_unique_forward.bedGraph.gz --pathCoverageReverse=${pathGenomeBrowser}/${sample}/coverage_unique_reverse.bedGraph.gz  --pathOutputGenes=${pathResults}/${sample}/CoverageGenes_IntronBlocks_${annot}.txt" >> ${pathScripts}/bsub_script

	if [ ${scheduler} = "slurm" ]; then
	    sbatch  ${pathScripts}/bsub_script
	else
	    chmod a+x ${pathScripts}/bsub_script
	    ${pathScripts}/bsub_script
	fi
    fi
fi



######################################################################################################
 
