#!/bin/bash
#SBATCH --wckey=A1297:SATURNE
#SBATCH -N 21
#SBATCH --time 10
srun hostname > $SPARK_HOME/conf/slaves
start-all.sh
R CMD BATCH corrigeMAP.r
###spark-submit --total-executor-cores 60 --executor-memory 5G germe_aleatoire.r 100
stop-all.sh
