#!/bin/bash

export sample=$1
export path=$2
export scheduler=$3
export ncores=$4

#############################################################################

export pathRNASeq=${path}/data/RNASeq
export pathDocs=${path}/docs
export pathResults=${path}/results/hisat2_alignments
export pathGenomeIndexes=${path}/results/hisat2_index
export pathScripts=${path}/scripts/RNASeq_alignments

export pathIndex=${pathGenomeIndexes}/genome_sequence

#############################################################################

if [ -e ${pathResults}/${sample}/accepted_hits.sam.gz ]||[ -e ${pathResults}/${sample}/accepted_hits.bam ]||[ -e ${pathResults}/${sample}/accepted_hits.sam ]||[ -e ${pathResults}/${sample}/accepted_hits_clean.sam ]||[ -e ${pathResults}/${sample}/accepted_hits_clean.sam.gz ]; then
    echo "already done (or ongoing)"
else
    
    export pathR1=""
    export pathR2=""
    
    export nbR1=0
    export nbR2=0
    export nbR=0
    
    for file in `grep $'\t'${sample}$ ${pathDocs}/RNASeq_Files.txt | cut -f 2 `
    do    
	if [ -e ${pathRNASeq}/${file}.fastq.gz ]; then
	    echo "already on /beegfs/"
	else
	    echo "downloading file from iRods"
	    iget MERIC/data/RNASeq/${file}.fastq.gz ${pathRNASeq}/
	fi
	
	export nbR=$((nbR+1))
    done
    
    for file in `grep $'\t'${sample}$ ${pathDocs}/RNASeq_Files.txt | cut -f 2 | grep R1`
    do
	if [ -e ${pathRNASeq}/${file}.fastq.gz ]; then
	    export pathR1=${pathRNASeq}/${file}.fastq.gz,${pathR1}
	    export nbR1=$((nbR1+1))
	fi
    done
    
    for file in `grep $'\t'${sample}$ ${pathDocs}/RNASeq_Files.txt | cut -f 2 | grep R2`
    do
	if [ -e ${pathRNASeq}/${file}.fastq.gz ]; then
	    export pathR2=${pathRNASeq}/${file}.fastq.gz,${pathR2}
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
    
    
    if [ -e ${pathResults}/${sample} ]; then
    	echo "dir output already there"
    else
    	mkdir ${pathResults}/${sample}
    fi
    
    export strand="--rna-strandness R"
    
    echo "#!/bin/bash" > ${pathScripts}/bsub_script_hisat

    if [ ${scheduler} = "slurm" ]; then
	echo "#SBATCH --job-name=hisat_${sample}" >>  ${pathScripts}/bsub_script_hisat
	echo "#SBATCH --output=${pathScripts}/std_out_${sample}" >>  ${pathScripts}/bsub_script_hisat
	echo "#SBATCH --error=${pathScripts}/std_err_${sample}" >>  ${pathScripts}/bsub_script_hisat
	echo "#SBATCH --cpus-per-task=${ncores}" >>  ${pathScripts}/bsub_script_hisat ## ${ncores} CPU
	echo "#SBATCH --time=24:00:00" >>  ${pathScripts}/bsub_script_hisat ## 24 hours
	echo "#SBATCH --mem=4G" >>  ${pathScripts}/bsub_script_hisat ## 5g per CPU
    fi	

    if [ ${nbR2} = 0 ]; then
	echo "single-end"
	echo "hisat2 --seed 19 -p ${ncores} -x ${pathIndex} -U ${pathR1} -S ${pathResults}/${sample}/accepted_hits.sam ${strand} --max-intronlen 1000000 --dta-cufflinks --no-unal --met-file ${pathResults}/${sample}/metrics.txt  --novel-splicesite-outfile ${pathResults}/${sample}/novel_splicesites.txt >& ${pathResults}/${sample}/align_summary.txt">> ${pathScripts}/bsub_script_hisat
    else
	if [ ${nbR2} = ${nbR1} ]; then
	    echo "paired-end"
	    echo "hisat2 --seed 19 -p ${ncores} -x ${pathIndex} -1 ${pathR1} -2 ${pathR2} -S ${pathResults}/${sample}/accepted_hits.sam ${strand} --max-intronlen 1000000 --dta-cufflinks --no-unal --met-file ${pathResults}/${sample}/metrics.txt  --novel-splicesite-outfile ${pathResults}/${sample}/novel_splicesites.txt >& ${pathResults}/${sample}/align_summary.txt">> ${pathScripts}/bsub_script_hisat
	fi
    fi
  
    for file in `grep $'\t'${sample}$ ${pathDocs}/RNASeq_Files.txt | cut -f 2 `
    do
	echo "rm ${pathRNASeq}/${file}.fastq.gz" >>  ${pathScripts}/bsub_script_hisat
    done

    echo "perl ${pathScripts}/cleanup.sam.pl --pathInput=${pathResults}/${sample}/accepted_hits.sam --pathOutput=${pathResults}/${sample}/accepted_hits_clean.sam " >> ${pathScripts}/bsub_script_hisat
    
    echo "rm ${pathResults}/${sample}/accepted_hits.sam " >>  ${pathScripts}/bsub_script_hisat
    
    echo "gzip ${pathResults}/${sample}/accepted_hits_clean.sam " >> ${pathScripts}/bsub_script_hisat
    
    if [ ${scheduler} = "slurm" ]; then
	sbatch ${pathScripts}/bsub_script_hisat
    else
	chmod +x ${pathScripts}/bsub_script_hisat
	${pathScripts}/bsub_script_hisat
    fi
    
fi

#############################################################################
