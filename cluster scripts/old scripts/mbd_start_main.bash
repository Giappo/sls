#!/bin/bash
#SBATCH --time=2:57:58

lambda_M=$1
mu_M=$2
lambda_S=$3
mu_S=$4
max_sims=$5

cd /home/$USER/mbd_like/

chmod +x install_packages.bash
./install_packages.bash --output=testinst.out
#sbatch install_packages.bash

sleep 60

rm testinst.out

if [[ -d sims/$lambda_M-$mu_M-$lambda_S-$mu_S/data ]]; then
  cd sims/$lambda_M-$mu_M-$lambda_S-$mu_S/
  sbatch /home/$USER/mbd_like/mbd_start_ML.bash $max_sims
else
  mkdir -p ./sims/$lambda_M-$mu_M-$lambda_S-$mu_S/data/
  sleep 1
  cd sims/$lambda_M-$mu_M-$lambda_S-$mu_S/
  
  echo "library(MBD)" > zzz_sim.R
  echo "args = as.numeric(commandArgs(TRUE))" >> zzz_sim.R
  echo "MBD:::mbd_sim_dataset(sim_pars=c(args[1],args[2],args[3],args[4]),max_sims=args[5],tips_interval=c(0,70),cond=1)" >> zzz_sim.R
  #module load R/3.3.1-foss-2016a
  module load R/3.4.4-foss-2018a-X11-20180131
  Rscript zzz_sim.R $lambda $mu $nu $q $max_sims
  
  wait zzz_sim
  
  sbatch /home/$USER/mbd_like/mbd_start_ML.bash $max_sims
fi


