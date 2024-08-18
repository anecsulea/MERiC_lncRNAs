#!/bin/bash

export sample=$1

############################################################################

export pathDocs=../../docs
export pathData=../../data/RNASeq
export pathResults=../../data_for_publication/EGA_submission/encrypted_files
export pathEGA=~/Tools/EGA-Cryptor-2.0.0
export iRods=MERIC/data/RNASeq/

mkdir -p ${pathResults}

############################################################################

for file in `grep $'\t'${sample}$ ${pathDocs}/RNASeq_Files.txt | cut -f 2`
do
    if [ -e ${pathData}/${file}.fastq.gz ]; then
	echo "ok, file is there"
    else
	cd ${pathData}
	iget ${iRods}/${file}.fastq.gz
	cd -
    fi
    

    if [ -e ${pathData}/encrypted_files/${file}.fastq.gz.gpg ]; then
	echo "already done"
    else
	java -jar ${pathEGA}/ega-cryptor-2.0.0.jar -t 4 -i ${pathData}/${file}.fastq.gz -o ${pathResults}/

	if [ -e ${pathData}/encrypted_files/${file}.fastq.gz.gpg ]; then
	    echo "ok"
	    rm ${pathData}/${file}.fastq.gz
	fi
    fi
done

############################################################################
	    
