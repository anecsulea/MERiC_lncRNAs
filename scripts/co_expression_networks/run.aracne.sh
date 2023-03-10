#!/bin/bash

export exp_dataset=$1
export TF_dataset=$2
export path=$3
export pathTools=$4
export threads=$5
export scheduler=$6
export seed=$7

####################################################################################

export pathResults=${path}/results/co_expression_networks/${exp_dataset}
export pathScripts=${path}/scripts/co_expression_networks

####################################################################################


echo "#!/bin/bash " > ${pathScripts}/bsub_script

if [ ${scheduler} = "slurm" ]; then
    echo "#SBATCH --job-name=aracne_${i}" >>  ${pathScripts}/bsub_script
    echo "#SBATCH --output=${pathScripts}/std_out_aracne_${seed}" >>  ${pathScripts}/bsub_script
    echo "#SBATCH --error=${pathScripts}/std_err_aracne_${seed}" >>  ${pathScripts}/bsub_script
    echo "#SBATCH --cpus-per-task=${threads}" >>  ${pathScripts}/bsub_script 
    echo "#SBATCH --time=4:00:00" >>  ${pathScripts}/bsub_script 
    echo "#SBATCH --mem=20G" >>  ${pathScripts}/bsub_script
fi

echo "java -Xmx5G -jar ${pathTools}/ARACNe-AP/dist/aracne.jar -e ${pathResults}/expression_data.txt -o ${pathResults}/${TF_dataset} --tfs ${pathResults}/${TF_dataset}.txt --pvalue 1E-8 --seed ${seed} --threads=${threads}"  >>  ${pathScripts}/bsub_script

if [ ${scheduler} = "slurm" ]; then
    sbatch ${pathScripts}/bsub_script
fi



####################################################################################
