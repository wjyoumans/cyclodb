#!/bin/bash
#SBATCH --job-name=cdb-3-50
#SBATCH --output=log/%j.out  # output and error messages go to this file
#SBATCH --mail-user=wyoumans@usf.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=rra
#SBATCH --qos=rra_guest
#SBATCH --ntasks=1
#SBATCH --array=3-3 #50
#SBATCH --mem-per-cpu=45GB
#SBATCH --time=07-00:00:00  # Wall Clock time (dd-hh:mm:ss) [max of 14 days]

# ntasks : number of parallel processes
# ntasks-per-node : number of processes to run on one node
# mem : specify the real memory required per node
# mem-per-cpu : minimum memory required per allocated CPU

# first: compute for all conductors up to 50

echo "$SLURM_JOB_ID $SLURM_JOB_NAME $SLURM_ARRAY_TASK_ID" >> monitor
time echo "cyclodata($SLURM_ARRAY_TASK_ID, prec=500, filename=\"$SLURM_ARRAY_TASK_ID\")" | gp -q -D path=.:abelianbnf -D parisize=45000000000 cyclodata.gp
