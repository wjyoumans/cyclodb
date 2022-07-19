#!/bin/bash
#SBATCH --job-name=cdb-real-101-200
#SBATCH --output=log/%j.out  # output and error messages go to this file
#SBATCH --mail-user=wyoumans@usf.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=rra
#SBATCH --qos=rra_guest
#SBATCH --ntasks=1
#SBATCH --array=101-200
#SBATCH --mem-per-cpu=185GB
#SBATCH --time=35-00:00:00  # Wall Clock time (dd-hh:mm:ss) [max of 14 days]

# ntasks : number of parallel processes
# ntasks-per-node : number of processes to run on one node
# mem : specify the real memory required per node
# mem-per-cpu : minimum memory required per allocated CPU

time echo "realcyclodata($SLURM_ARRAY_TASK_ID, prec=500, filename=\"out/$SLURM_ARRAY_TASK_ID\")" | gp -q -D path=.:abelianbnf -D parisize=185000000000 cyclodata.gp
