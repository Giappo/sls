#!/bin/bash
#SBATCH --time=00:57:58 --partition=gelifes

s = 1
lambda_M=$1
mu_M=$2
lambda_S=$3
mu_S=$4
cond=$5
max_sims=$6
model=$7
function="sls::loglik_"$model

echo "library(sls)" > test.R
echo "args0 = as.numeric(commandArgs(TRUE))" >> test.R
echo "args = commandArgs(TRUE)" >> test.R

#echo "args0 = as.numeric(commandArgs(TRUE))" >> test.R
#echo "args <- as.numeric(args0)" >> test.R
#echo "args[is.na(args)] <- args0[is.na(args)]" >> test.R
echo "for (i in seq_along(args0)) {print(args0[i])}" >> test.R
echo "for (i in seq_along(args)) {print(args[i])}" >> test.R

Rscript test.R $s $lambda_M $mu_M $lambda_S $mu_S $cond $function