#!/bin/bash
#SBATCH --time=7:57:58

Nsims=$1

cd /home/$USER/sls/

sbatch install_packages.bash --output=testinst.out
wait install_packages.bash #sleep 30
rm testinst.out

echo "library(sls)" > sls_determine_length_datasets.R
echo "load_all_data(the.environment = environment()); data.sets <- ls(pattern = 'dataset_',envir = environment())" >> sls_determine_length_datasets.R
echo "L <- length(data.sets)" >> sls_determine_length_datasets.R
echo "sink('LengthDataset.txt')" >> sls_determine_length_datasets.R
echo "cat(paste0('LENGTH=',L))" >> sls_determine_length_datasets.R
echo "sink()" >> sls_determine_length_datasets.R
echo "q(save = 'no')" >> sls_determine_length_datasets.R
module load R/3.3.1-foss-2016a
Rscript sls_determine_length_datasets.R

sleep 5

source "LengthDataset.txt"
L=$LENGTH
rm "LengthDataset.txt"
rm "sls_determine_length_datasets.R"

sleep 1

rm runtest.R
#echo "library(expoRkit)" > zzz_ML.R
echo "library(sls)" >> runtest.R
echo "args = as.numeric(commandArgs(TRUE))" >> runtest.R
echo "sls:::test_likelihood_formula2(s=args[1], Nsims=args[2], lik_function = lik_custom, sim_function = sim_custom4)" >> runtest.R

for((s = 1; s <= L; s++))
do

#echo "#!/bin/bash" > sls_job$s
echo '#!/bin/bash' > sls_job$s
echo "#SBATCH --time=71:59:00" >> sls_job$s
echo "module load R/3.3.1-foss-2016a" >> sls_job$s
echo "Rscript runtest.R $s $Nsims" >> sls_job$s
echo "rm sls_job$s" >> sls_job$s

sbatch --partition=regular --mem=9GB --job-name=sls$s --mail-type=FAIL,TIME_LIMIT --mail-user=glaudanno@gmail.com sls_job$s

done








#echo "library(sls)" > zzz_sim.R
#echo "args = as.numeric(commandArgs(TRUE))" >> zzz_sim.R
#echo "MBD:::mbd_sim_dataset2(sim_pars=c(args[1],args[2],args[3],args[4]),max_sims=args[5],tips_interval=c(0,70),cond=1)" >> zzz_sim.R
#module load R/3.3.1-foss-2016a
#Rscript sls_determine_length_datasets.R 
#Rscript zzz_sim.R $lambda $mu $nu $q $max_sims
  
#wait zzz_sim
 
#sbatch /home/$USER/sls/sls_start_test.bash $Nsims