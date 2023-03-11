#!/bin/bash

export sample=$1
export path=$2
export scheduler=$3
export ncores=$4

#############################################################################

export pathResults=${path}/results/hisat2_alignments
export pathScripts=${path}/scripts/RNASeq_alignments

#############################################################################

if [ -e ${pathResults}/${sample}/accepted_hits.bam ]; then
    echo "already done"
    exit
fi

#############################################################################

echo "#!/bin/bash" > ${pathScripts}/bsub_script_sort

if [ ${scheduler} = "slurm" ]; then
    echo "#SBATCH --job-name=sort_${sample}" >>  ${pathScripts}/bsub_script_hisat
    echo "#SBATCH --output=${pathScripts}/std_out_sort_${sample}" >>  ${pathScripts}/bsub_script_hisat
    echo "#SBATCH --error=${pathScripts}/std_err_sort_${sample}" >>  ${pathScripts}/bsub_script_hisat
    echo "#SBATCH --cpus-per-task=${ncores}" >>  ${pathScripts}/bsub_script_hisat ## ${ncores} CPU
    echo "#SBATCH --time=24:00:00" >>  ${pathScripts}/bsub_script_hisat ## 24 hours
    echo "#SBATCH --mem=4G" >>  ${pathScripts}/bsub_script_hisat ## 5g per CPU
fi	

echo "samtools sort -m 5G -@ ${ncores} -o ${pathResults}/${sample}/accepted_hits.bam -O bam ${pathResults}/${sample}/accepted_hits.sam" >> ${pathScripts}/bsub_script_sort
echo "if [ -e  ${pathResults}/${sample}/accepted_hits.bam ]; then" >> ${pathScripts}/bsub_script_sort  
echo "rm ${pathResults}/${sample}/accepted_hits.sam" >> ${pathScripts}/bsub_script_sort  
echo "fi"  >> ${pathScripts}/bsub_script_sort  

#############################################################################

if [ ${scheduler} = "slurm" ]; then
    sbatch ${pathScripts}/bsub_script_sort  
else
    chmod a+x ${pathScripts}/bsub_script_sort
    ${pathScripts}/bsub_script_sort
fi

#############################################################################
