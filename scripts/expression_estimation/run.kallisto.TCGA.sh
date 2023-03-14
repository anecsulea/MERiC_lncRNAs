#!/bin/bash

export sample=$1
export annot=$2
export path=$3
export pathTools=$4
export nthreads=$5

####################################################################################

export pathDocs=${path}/data/TCGA
export pathRNASeq=${path}/data/TCGA/fastq_files
export pathResults=${path}/results/expression_estimation_TCGA
export pathIndexes=${path}/results/kallisto_indexes
export pathScripts=${path}/scripts/expression_estimation

####################################################################################

if [ -e ${pathResults}/${annot}/${sample}/abundance.tsv ]; then
    echo "already done"
else

    export pathR1=""
    export pathR2=""

    export nbR1=0
    export nbR2=0
    export nbR=0

    export dir=`grep $'\t'${sample}$ ${pathDocs}/lihc_rnaseq_sample_info.txt | cut -f 1 `
    export file=`grep $'\t'${sample}$ ${pathDocs}/lihc_rnaseq_sample_info.txt | cut -f 2 `

    export pathR1=${pathRNASeq}/${dir}/${file}_R1.fq.gz
    export pathR2=${pathRNASeq}/${dir}/${file}_R2.fq.gz

    if [ -e ${pathR1} ]; then
	echo "already in ${pathRNASeq}"
    else
	echo "downloading file from iRods"
	iget MERIC/data/TCGA/fastq_files/${dir}/${file}_R1.fq.gz ${pathRNASeq}/${dir}/
    fi
    
     if [ -e ${pathR2} ]; then
	 echo "already in ${pathRNASeq}"
     else
	echo "downloading file from iRods"
	iget MERIC/data/TCGA/fastq_files/${dir}/${file}_R2.fq.gz ${pathRNASeq}/${dir}/
     fi

    #############################################################################

    if [ -e ${pathResults}/${annot}/${sample} ]; then
	echo "path exists"
    else
	mkdir -p ${pathResults}/${annot}/${sample}
    fi

     #############################################################################
    
    echo "#!/bin/bash" >  ${pathScripts}/bsub_script
        
    echo "#SBATCH --job-name=kallisto_${sample}" >>  ${pathScripts}/bsub_script
    echo "#SBATCH --partition=normal" >>  ${pathScripts}/bsub_script
    echo "#SBATCH --output=${pathScripts}/std_out_kallisto_${sample}" >>  ${pathScripts}/bsub_script
    echo "#SBATCH --error=${pathScripts}/std_err_kallisto_${sample}" >>  ${pathScripts}/bsub_script
    echo "#SBATCH --cpus-per-task=${nthreads}" >>  ${pathScripts}/bsub_script ## ${ncores} CPU
    echo "#SBATCH --time=2:00:00" >>  ${pathScripts}/bsub_script ## 8 hours
    echo "#SBATCH --mem=25G" >>  ${pathScripts}/bsub_script ## 15g per CPU
    
    echo "singularity exec -B ${path} ${pathTools}/kallisto.sif kallisto quant  --bias -t ${nthreads} -o ${pathResults}/${annot}/${sample}/ --index ${pathIndexes}/${annot} ${pathR1} ${pathR2}" >> ${pathScripts}/bsub_script

    ## cleanup if everything went well
    echo "if [ -e ${pathResults}/${annot}/${sample}/abundance.tsv ]; then" >>  ${pathScripts}/bsub_script
    echo "rm ${pathR1}" >>  ${pathScripts}/bsub_script
    echo "rm ${pathR2}" >>  ${pathScripts}/bsub_script
    echo "fi" >>  ${pathScripts}/bsub_script

    sbatch ${pathScripts}/bsub_script
fi

####################################################################################
