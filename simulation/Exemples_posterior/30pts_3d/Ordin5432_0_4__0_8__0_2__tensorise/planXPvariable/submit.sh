#!/bin/bash
#SBATCH --wckey=P11E1:CARMEL_3D
#SBATCH -N 21
#SBATCH --time 60
srun hostname > $SPARK_HOME/conf/slaves
start-all.sh
R CMD BATCH germe_aleatoire.r
###spark-submit --total-executor-cores 60 --executor-memory 5G germe_aleatoire.r 100
stop-all.sh
