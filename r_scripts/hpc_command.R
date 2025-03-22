#Randi script 
#Log in 
ssh beiw@randi.cri.uchicago.edu

#remove 
rm /gpfs/data/chiu-lab/bw_folder/dissertation/5hmc_tx/interation_100_allrace_newMM/res/final_*
#remove logs
rm /gpfs/data/chiu-lab/bw_folder/dissertation/5hmc_tx/interation_100_allrace_newMM/logs/res_*
  
#submit jobs 
for i in `seq 1 100`; do sbatch --export=ARGS1=${i} /gpfs/data/chiu-lab/bw_folder/dissertation/5hmc_tx/interation_100_allrace_newMM/sbatch/feature_selection.sh; done

#check status 
squeue -u beiw 
#check logs 
less /gpfs/data/chiu-lab/bw_folder/dissertation/5hmc_tx/interation_100_allrace_newMM/logs/res_*

##Direct to chiu lab: 
cd /gpfs/data/chiu-lab

##Chiu-lab path
path=/gpfs/data/chiu-lab/
  