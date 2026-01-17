#!/bin/bash -l

#SBATCH --output=/home/projects/22126_NGS/projects/group15/logs/%x_o.%A
#SBATCH --error=/home/projects/22126_NGS/projects/group15/logs/%x_e.%A
#SBATCH -c 4
#SBATCH --mem=10G
#SBATCH --time=02:00:00 

PROJECT_HOME_DIR=${PROJECT_HOME_DIR:-/home/projects/22126_NGS/projects/group15}
PROJECT_SCRIPTS_DIR=$PROJECT_HOME_DIR/scripts

preprocess=$PROJECT_SCRIPTS_DIR/preprocess.sh

cd $PROJECT_HOME_DIR
echo "Running $preprocess with params: $@"

$preprocess $@
