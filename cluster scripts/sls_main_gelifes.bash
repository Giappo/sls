#!/bin/bash
#SBATCH --time=00:57:58 --partition=gelifes

lambda_M=$1
mu_M=$2
lambda_S=$3
mu_S=$4
cond=$5
max_sims=$6
model=$7
function="sls::loglik_"$model

cd /home/$USER/sls/

chmod +x install_packages.bash
./install_packages.bash --output=testinst.out
#sbatch install_packages.bash

sleep 40

rm testinst.out

if [[ -d sims/$lambda_M-$mu_M-$lambda_S-$mu_S-$cond/data ]]; then
#if [[ -d sims/$model1-vs-$model2/$lambda_M-$mu_M-$lambda_S-$mu_S-$cond/data ]]; then
  cd sims/$lambda_M-$mu_M-$lambda_S-$mu_S-$cond/
else
  mkdir -p ./sims/$lambda_M-$mu_M-$lambda_S-$mu_S-$cond/data/
  #mkdir -p ./sims/$model1-vs-$model2/$lambda_M-$mu_M-$lambda_S-$mu_S-$cond/data/
  sleep 1
  cd sims/$lambda_M-$mu_M-$lambda_S-$mu_S-$cond/
fi

echo "rm -rfv errors/*"
#sbatch /home/$USER/mbd_like/install_packages.bash & sleep 60

echo "library(sls)" > zzz_ML_gelifes.R
echo "args = as.numeric(commandArgs(TRUE))" >> zzz_ML_gelifes.R
#echo "sls:::sls_ML_cluster(s=args[1],simpars=c(args[2],args[3],args[4],args[5]),cond=args[6],t_d=4,fun=eval(parse(text = paste0('sls::loglik_', toString(args[7])))))" >> zzz_ML_gelifes.R
#echo "sls:::sls_ML_cluster(s=args[1],simpars=c(args[2],args[3],args[4],args[5]),cond=args[6],t_d=4,fun=args[7])" >> zzz_ML_gelifes.R
#echo "sls:::sls_ML_cluster(s=args[1],simpars=c(args[2],args[3],args[4],args[5]),cond=args[6],t_d=4,fun=eval(parse(text = args[7])))" >> zzz_ML_gelifes.R
echo "sls:::sls_ML_cluster(s=args[1],simpars=c(args[2],args[3],args[4],args[5]),cond=args[6],t_d=4,fun=eval(args[7]))" >> zzz_ML_gelifes.R

for((s = 1; s <= max_sims; s++))
do

echo "#!/bin/bash" > SLSjob$s$lambda_M$mu_M$lambda_S$mu_S$cond$model
echo "#SBATCH --time=169:59:00" >> SLSjob$s$lambda_M$mu_M$lambda_S$mu_S$cond$model
echo "module load R/3.4.4-foss-2018a-X11-20180131" >> SLSjob$s$lambda_M$mu_M$lambda_S$mu_S$cond$model
#echo "Rscript zzz_ML_gelifes.R $s $lambda_M $mu_M $lambda_S $mu_S $cond $model" >> SLSjob$s$lambda_M$mu_M$lambda_S$mu_S$cond$model
echo "Rscript zzz_ML_gelifes.R $s $lambda_M $mu_M $lambda_S $mu_S $cond $function" >> SLSjob$s$lambda_M$mu_M$lambda_S$mu_S$cond$model
echo "rm SLSjob$s$lambda_M$mu_M$lambda_S$mu_S$cond$model" >> SLSjob$s$lambda_M$mu_M$lambda_S$mu_S$cond$model

#NEVER ASK FOR MORE THAN 9GB OF MEMORY!
#sbatch --partition=regular --mem=9GB --job-name=ML$s$lambda_M$mu_M$lambda_S$mu_S$cond$model --mail-type=FAIL,TIME_LIMIT --mail-user=glaudanno@gmail.com SLSjob$s$lambda_M$mu_M$lambda_S$mu_S$cond$model
sbatch --partition=gelifes --mem=9GB --job-name=ML$s$lambda_M$mu_M$lambda_S$mu_S$cond$model --mail-type=FAIL,TIME_LIMIT --mail-user=glaudanno@gmail.com SLSjob$s$lambda_M$mu_M$lambda_S$mu_S$cond$model

done

