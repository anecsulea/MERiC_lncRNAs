#!/bin/bash

export sample=$1
export annot=$2
export path=$3
export pathTools=$4
export nthreads=$5

####################################################################################

export pathDocs=${path}/docs
export pathRNASeq=${path}/data/RNASeq
export pathResults=${path}/results/expression_estimation
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

    for file in `grep $'\t'${sample}$ ${pathDocs}/RNASeq_Files.txt | cut -f 2 `
    do
	if [ -e ${pathRNASeq}/${file}.fastq.gz ]; then
	    echo "already in ${pathRNASeq}"
	else
	    echo "downloading file from iRods"
	    iget MERIC/data/RNASeq/${file}.fastq.gz ${pathRNASeq}/
	fi

	export nbR=$((nbR+1))
    done

    for file in `grep $'\t'${sample}$ ${pathDocs}/RNASeq_Files.txt | cut -f 2 | grep R1`
    do
	if [ -e ${pathRNASeq}/${file}.fastq.gz ]; then
	    export pathR1=${pathRNASeq}/${file}.fastq.gz" "${pathR1}
	    export nbR1=$((nbR1+1))
	fi
    done

    for file in `grep $'\t'${sample}$ ${pathDocs}/RNASeq_Files.txt | cut -f 2 | grep R2`
    do
	if [ -e ${pathRNASeq}/${file}.fastq.gz ]; then
	    export pathR2=${pathRNASeq}/${file}.fastq.gz" "${pathR2}
	    export nbR2=$((nbR2+1))
	fi
    done

    export nbth=`grep -c $'\t'${sample}$ ${pathDocs}/RNASeq_Files.txt`

    if [ ${nbth} = ${nbR} ]&&[ ${nbR} -gt 0 ]; then
	echo "all files found for "${sample} ${pathR1} ${pathR2} ${nbR1} ${nbR2}
    else
	echo "not doing anything for "${sample}", files are not there"
	echo ${nbth} ${nbR} ${nbR1} ${nbR2}
	exit
    fi

    #############################################################################

    if [ -e ${pathResults}/${annot}/${sample} ]; then
	echo "path exists"
    else
	mkdir -p ${pathResults}/${annot}/${sample}
    fi

    
    echo "#!/bin/bash " > ${pathScripts}/bsub_script
    echo "#SBATCH --job-name=kallisto_${sample}" >>  ${pathScripts}/bsub_script
    echo "#SBATCH --partition=normal" >>  ${pathScripts}/bsub_script
    echo "#SBATCH --output=${pathScripts}/std_out_kallisto_${sample}" >>  ${pathScripts}/bsub_script
    echo "#SBATCH --error=${pathScripts}/std_err_kallisto_${sample}" >>  ${pathScripts}/bsub_script
    echo "#SBATCH --cpus-per-task=${nthreads}" >>  ${pathScripts}/bsub_script ## ${ncores} CPU
    echo "#SBATCH --time=1:00:00" >>  ${pathScripts}/bsub_script ## 8 hours
    echo "#SBATCH --mem=15G" >>  ${pathScripts}/bsub_script ## 15g per CPU
    
    echo "singularity exec -B ${path} ${pathTools}/kallisto.sif kallisto quant --single -l 200.0 -s 25 --bias -t ${nthreads} --rf-stranded -o ${pathResults}/${annot}/${sample} --index ${pathIndexes}/${annot} ${pathR1} " >> ${pathScripts}/bsub_script

    echo "if [ -e ${pathResults}/${annot}/${sample}/abundance.tsv ]; then" >>  ${pathScripts}/bsub_script

    for file in `grep $'\t'${sample}$ ${pathDocs}/RNASeq_Files.txt | cut -f 2 `
    do
	echo "rm ${pathRNASeq}/${file}.fastq.gz" >>  ${pathScripts}/bsub_script
    done

    echo "fi" >>  ${pathScripts}/bsub_script

    sbatch ${pathScripts}/bsub_script
   
fi

####################################################################################
