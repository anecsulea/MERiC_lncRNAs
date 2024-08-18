#!/bin/bash

export sample=$1
export path=$2
export pathTools=$3

############################################################################

export pathScripts=${path}/scripts/prepare_EGA_submission
export pathDocs=${path}/docs
export pathData=${path}/data/RNASeq
export pathResults=${path}/data_for_publication/EGA_submission/encrypted_files
export pathEGA=${pathTools}/EGA-Cryptor-2.0.0

export iRods=MERIC/data/RNASeq/

############################################################################

mkdir -p ${pathResults}

############################################################################

for file in `grep $'\t'${sample}$ ${pathDocs}/RNASeq_Files.txt | cut -f 2`
do
    if [ -e ${pathData}/encrypted_files/${file}.fastq.gz.gpg ]; then
	echo "already done"
    else
	if [ -e ${pathData}/${file}.fastq.gz ]; then
	    echo "ok, file is there"
	else
	    cd ${pathData}
	    iget ${iRods}/${file}.fastq.gz
	    cd ${pathScripts}
	fi

	echo "#!/bin/bash " > ${pathScripts}/bsub_script
	echo "#SBATCH --job-name=encrypt_${file}" >>  ${pathScripts}/bsub_script
	echo "#SBATCH --partition=normal" >>  ${pathScripts}/bsub_script
	echo "#SBATCH --output=${pathScripts}/std_out_encrypt_${file}" >>  ${pathScripts}/bsub_script
	echo "#SBATCH --error=${pathScripts}/std_err_encrypt_${file}" >>  ${pathScripts}/bsub_script
	echo "#SBATCH --cpus-per-task=4" >>  ${pathScripts}/bsub_script ## 4 CPU
	echo "#SBATCH --time=1:00:00" >>  ${pathScripts}/bsub_script ## 1 hour
	echo "#SBATCH --mem=10G" >>  ${pathScripts}/bsub_script ## per CPU
    
	echo "java -jar ${pathEGA}/ega-cryptor-2.0.0.jar -t 4 -i ${pathData}/${file}.fastq.gz -o ${pathResults}/" >>  ${pathScripts}/bsub_script

	echo "if [ -e ${pathResults}/${file}.fastq.gz.gpg ]; then" >>  ${pathScripts}/bsub_script
	echo "echo \"ok\"" >>  ${pathScripts}/bsub_script
	echo "rm ${pathData}/${file}.fastq.gz" >>  ${pathScripts}/bsub_script
	echo "fi" >>  ${pathScripts}/bsub_script

	sbatch ${pathScripts}/bsub_script
    fi
done

############################################################################
	    
