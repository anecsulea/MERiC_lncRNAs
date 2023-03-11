#!/bin/bash

export sample=$1
export path=$2
export scheduler=$3

###################################################################################

export pathAlignments=${path}/results/hisat2_alignments
export pathScripts=${path}/scripts/genome_browser

###################################################################################

if [ -e ${pathAlignments}/${sample}/accepted_hits.bam ]; then
    
    if [ -e ${pathAlignments}/${sample}/accepted_hits_unique_forward.bam ]&&[ -e ${pathAlignments}/${sample}/accepted_hits_unique_reverse.bam ]; then
	echo "already done"
    else
	echo ${sample}
	
	echo "#!/bin/bash" >  ${pathScripts}/bsub_script_unique_${sample}
	
	if [ ${scheduler} = "slurm" ]; then
	    echo "#SBATCH --job-name=hisat_${sample}" >>  ${pathScripts}/bsub_script_unique_${sample}
	    echo "#SBATCH --output=${pathScripts}/std_out_${sample}" >>  ${pathScripts}/bsub_script_unique_${sample}
	    echo "#SBATCH --error=${pathScripts}/std_err_${sample}" >>  ${pathScripts}/bsub_script_unique_${sample}
	    echo "#SBATCH --cpus-per-task=1" >> ${pathScripts}/bsub_script_unique_${sample}
	    echo "#SBATCH --time=24:00:00" >>  ${pathScripts}/bsub_script_unique_${sample}
	    echo "#SBATCH --mem=4G" >>  ${pathScripts}/bsub_script_unique_${sample}
	fi
	
	## unique reads
	
	echo "perl ${pathScripts}/extract.unique.reads.pl --pathAllReads=${pathAlignments}/${sample}/accepted_hits.bam --pathUniqueReads=${pathAlignments}/${sample}/accepted_hits_unique.sam " >>  ${pathScripts}/bsub_script_unique_${sample}
	
	echo "grep ^@ ${pathAlignments}/${sample}/accepted_hits_unique.sam > ${pathAlignments}/${sample}/accepted_hits_unique_forward.sam" >>  ${pathScripts}/bsub_script_unique_${sample}
	echo "cp ${pathAlignments}/${sample}/accepted_hits_unique_forward.sam ${pathAlignments}/${sample}/accepted_hits_unique_reverse.sam" >>  ${pathScripts}/bsub_script_unique_${sample} 
	
	
	## forward strand
	
	echo "grep $'\t'XS:A:+ ${pathAlignments}/${sample}/accepted_hits_unique.sam >> ${pathAlignments}/${sample}/accepted_hits_unique_forward.sam" >>  ${pathScripts}/bsub_script_unique_${sample}
	echo "samtools view -Sb  -o ${pathAlignments}/${sample}/accepted_hits_unique_forward.bam ${pathAlignments}/${sample}/accepted_hits_unique_forward.sam" >>  ${pathScripts}/bsub_script_unique_${sample}
	echo "rm ${pathAlignments}/${sample}/acc*forward*sam" >>  ${pathScripts}/bsub_script_unique_${sample}
	
	
	## reverse strand 
	
	echo "grep $'\t'XS:A:- ${pathAlignments}/${sample}/accepted_hits_unique.sam >> ${pathAlignments}/${sample}/accepted_hits_unique_reverse.sam" >>  ${pathScripts}/bsub_script_unique_${sample}
	
	echo "rm ${pathAlignments}/${sample}/accepted_hits_unique.sam " >>  ${pathScripts}/bsub_script_unique_${sample}
	
	echo "samtools view -Sb  -o ${pathAlignments}/${sample}/accepted_hits_unique_reverse.bam ${pathAlignments}/${sample}/accepted_hits_unique_reverse.sam" >>  ${pathScripts}/bsub_script_unique_${sample}
	echo "rm ${pathAlignments}/${sample}/acc*reverse*sam" >>  ${pathScripts}/bsub_script_unique_${sample}
	
	
	## cleanup 
	echo "rm ${pathAlignments}/${sample}/accepted_hits_unique.sam" >>  ${pathScripts}/bsub_script_unique_${sample}
	
	if [ ${scheduler} = "slurm" ]; then
	    sbatch ${pathScripts}/bsub_script_unique_${sample}
	else
	    chmod a+x ${pathScripts}/bsub_script_unique_${sample}
	    ${pathScripts}/bsub_script_unique_${sample}
	fi
    fi
fi


###################################################################################
