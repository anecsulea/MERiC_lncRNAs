#!/bin/bash

export sample=$1
export path=$2
export scheduler=$3

###################################################################################

export pathAlignments=${path}/results/hisat2_alignments
export pathResults=${path}/results/genome_browser
export pathScripts=${path}/scripts/genome_browser

###################################################################################

if [ -e ${pathAlignments}/${sample}/accepted_hits_unique_forward.bam ]&&[ -e ${pathAlignments}/${sample}/accepted_hits_unique_reverse.bam ]; then
    
    if [ -e ${pathResults}/${sample} ]; then
	echo "dir output already there"
    else
	mkdir ${pathResults}/${sample} 
    fi
    
    if [ -e ${pathResults}/${sample}/coverage_unique_reverse.bedGraph.gz ]&&[ -e ${pathResults}/${sample}/coverage_unique_forward.bedGraph.gz ]; then
	    echo "already done"
    else
	if [ -e ${pathResults}/${sample}/coverage_unique_reverse.bedGraph ]&&[ -e ${pathResults}/${sample}/coverage_unique_forward.bedGraph ]; then
	    echo "already done"
	else
	    echo ${sample}
	    
	    echo "#!/bin/bash" >  ${pathScripts}/bsub_script_coverage_${sample}
	    
	    if [ ${scheduler} = "slurm" ]; then
		echo "#SBATCH --job-name=coverage_${sample}" >>  ${pathScripts}/bsub_script_coverage_${sample}
		echo "#SBATCH --output=${pathScripts}/std_out_coverage_${sample}" >>  ${pathScripts}/bsub_script_coverage_${sample}
		echo "#SBATCH --error=${pathScripts}/std_err_coverage_${sample}" >>  ${pathScripts}/bsub_script_coverage_${sample}
		echo "#SBATCH --cpus-per-task=1" >> ${pathScripts}/bsub_script_coverage_${sample}
		echo "#SBATCH --time=24:00:00" >>  ${pathScripts}/bsub_script_coverage_${sample}
		echo "#SBATCH --mem=12G" >>  ${pathScripts}/bsub_script_coverage_${sample}
	    fi
	    
	    ## forward strand
	    
	    echo "bedtools genomecov -bg -split -ibam ${pathAlignments}/${sample}/accepted_hits_unique_forward.bam > ${pathResults}/${sample}/coverage_unique_forward.bedGraph" >>  ${pathScripts}/bsub_script_coverage_${sample} 
	    
	    echo "gzip ${pathResults}/${sample}/cov*forward*" >>  ${pathScripts}/bsub_script_coverage_${sample}
	    
	    ## reverse strand 
	    
	    echo "bedtools genomecov -bg -split -ibam ${pathAlignments}/${sample}/accepted_hits_unique_reverse.bam > ${pathResults}/${sample}/coverage_unique_reverse.bedGraph " >>  ${pathScripts}/bsub_script_coverage_${sample}
	    
	    echo "gzip ${pathResults}/${sample}/cov*reverse*" >>  ${pathScripts}/bsub_script_coverage_${sample}
	    
	    if [ ${scheduler} = "slurm" ]; then
		sbatch ${pathScripts}/bsub_script_coverage_${sample}
	    else
	    	chmod a+x ${pathScripts}/bsub_script_coverage_${sample}
		${pathScripts}/bsub_script_coverage_${sample} 
	    fi
	fi
    fi
fi


###################################################################################
