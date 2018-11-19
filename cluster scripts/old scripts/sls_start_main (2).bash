#!/bin/bash
#SBATCH --time=00:57:58

lambda_M=$1
mu_M=$2
lambda_S=$3
mu_S=$4
cond=$5
max_sims=$6
#change the names of the job if you change parameters

cd /home/$USER/sls/

chmod +x install_packages.bash
./install_packages.bash --output=testinst.out
#sbatch install_packages.bash

sleep 40

rm testinst.out

if [[ -d sims/$lambda_M-$mu_M-$lambda_S-$mu_S-$cond/data ]]; then
  cd sims/$lambda_M-$mu_M-$lambda_S-$mu_S-$cond/
else
  mkdir -p ./sims/$lambda_M-$mu_M-$lambda_S-$mu_S-$cond/data/
  sleep 1
  cd sims/$lambda_M-$mu_M-$lambda_S-$mu_S-$cond/
fi

echo "rm -rfv errors/*"
#sbatch /home/$USER/mbd_like/install_packages.bash & sleep 60

echo "library(sls)" > zzz_ML.R
echo "args = as.numeric(commandArgs(TRUE))" >> zzz_ML.R
echo "sls:::sls_ML_cluster(s=args[1],simpars=c(args[2],args[3],args[4],args[5]),cond=args[6],t_d=4)" >> zzz_ML.R

for((s = 1; s <= max_sims; s++))
do

echo "#!/bin/bash" > SLSjob$s$lambda_M$mu_M$lambda_S$mu_S$cond
echo "#SBATCH --time=169:59:00" >> SLSjob$s$lambda_M$mu_M$lambda_S$mu_S$cond
#echo "module load R/3.3.1-foss-2016a" >> SLSjob$s
echo "module load R/3.4.4-foss-2018a-X11-20180131" >> SLSjob$s$lambda_M$mu_M$lambda_S$mu_S$cond
echo "Rscript zzz_ML.R $s $lambda_M $mu_M $lambda_S $mu_S $cond" >> SLSjob$s$lambda_M$mu_M$lambda_S$mu_S$cond
echo "rm SLSjob$s$lambda_M$mu_M$lambda_S$mu_S$cond" >> SLSjob$s$lambda_M$mu_M$lambda_S$mu_S$cond

#NEVER ASK FOR MORE THAN 9GB OF MEMORY!
sbatch --partition=regular --mem=9GB --job-name=ML$s$lambda_M$mu_M$lambda_S$mu_S$cond --mail-type=FAIL,TIME_LIMIT --mail-user=glaudanno@gmail.com SLSjob$s$lambda_M$mu_M$lambda_S$mu_S$cond

done

