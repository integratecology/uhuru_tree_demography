#!/bin/bash
#SBATCH --job-name=uhuru_trees     # name of the job
#SBATCH --partition=defq           # partition to be used (defq, gpu or intel)
#SBATCH --time=96:00:00            # walltime (up to 96 hours)
#SBATCH --nodes=1                  # number of nodes
#SBATCH --ntasks-per-node=1        # number of tasks (i.e. parallel processes) to be started
#SBATCH --cpus-per-task=1          # number of cpus required to run the script
#SBATCH --mem-per-cpu=32G	   # memory required for process
#SBATCH --array=0-2%10    	   # set number of total simulations and number that can run simultaneously	  


module load gcc
module load R

cd /home/alston92/proj/uhuru_tree_demography   # where executable and data is located

date
echo "Initiating script"


if [ -f results/stochastic_lambdas.csv ]; then
	echo "Results file already exists! continuing..."
else
	echo "creating results file for stochastic lambdas"
	echo "parameter,species,treatment,value" > results/stochastic_lambdas.csv
fi

if [ -f results/sensitivities.csv ]; then
        echo "Results file already exists! continuing..."
else
        echo "creating results file for stochastic lambdas"
        echo "parameter,species,treatment,vital_rate,value" > results/stochastic_lambdas.csv
fi

if [ -f results/elasticities.csv ]; then
        echo "Results file already exists! continuing..."
else
        echo "creating results file for stochastic lambdas"
        echo "parameter,species,treatment,vital_rate,value" > results/stochastic_lambdas.csv
fi


Rscript uhuru_tree_demography_analysis_full_script.R ${SLURM_ARRAY_TASK_ID}     # name of script
echo "Script complete"
date
