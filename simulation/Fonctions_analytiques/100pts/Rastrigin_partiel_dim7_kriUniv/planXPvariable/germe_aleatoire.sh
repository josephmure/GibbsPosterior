#!/bin/bash
#SBATCH --wckey=P10M2:CARMEL_3D
#SBATCH -N 21
#SBATCH --time 600
export SPARK_MASTER_HOST=atfront1.athos.hpc.edf.fr
export SPARK_MASTER_PORT=7077
export SPARK_CONF_DIR="`pwd`"
echo $SPARK_MASTER_HOST
echo $SPARK_CONF_DIR
srun hostname > $SPARK_CONF_DIR/slaves
start-all.sh
echo $SPARK_MASTER_HOST
echo $SPARK_CONF_DIR
R CMD BATCH germe_aleatoire.r
stop-all.sh
rm $SPARK_CONF_DIR/slaves
#start-master.sh
# ATTENTION, J'AI AJOUTE bash -i POUR FORCER L'EXECTUTION D'UN SHELL INTERACTIF : SANS CELA, LE SHELL EST EXECUTE EN MODE NON INTERACTIF ET NE LIT DONC PAS LE FICHIER .bashrc
#"${SPARK_HOME}/sbin/slaves.sh" cd "${SPARK_HOME}" \; bash -i "${SPARK_HOME}/sbin/start-slave.sh" "spark://$SPARK_MASTER_HOST:$SPARK_MASTER_PORT"
#R CMD BATCH germe_aleatoire.r
#"${SPARK_HOME}/sbin/slaves.sh" cd "${SPARK_HOME}" \; "${SPARK_HOME}/sbin"/stop-slave.sh
#stop-master.sh
