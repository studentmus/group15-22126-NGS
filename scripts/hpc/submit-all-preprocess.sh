#!/bin/bash

set -e
set -u
set -o pipefail

PROJECT_HOME_DIR=${PROJECT_HOME_DIR:-/home/projects/22126_NGS/projects/group15}
PROJECT_DATA_DIR=$PROJECT_HOME_DIR/data
LOG_DIR=$PROJECT_HOME_DIR/logs

mkdir -p -m 775 $LOG_DIR

log_file="$LOG_DIR/$(basename $0).log"

log() {
    local msg="[$(date +"%Y-%m-%d %H:%M:%S")] [$(basename $0)] $*" 
    echo "$msg" >> "$log_file" 2>&1
    echo "$msg"
}

srr_ids=$(cat $PROJECT_DATA_DIR/samples.csv | cut -d ',' -f1)

for srr_id in $srr_ids; do
    fastq_file=$PROJECT_DATA_DIR/raw/$srr_id.fastq.gz
    if [ -f $fastq_file ]; then
        log "Submitting job for $fastq_file"
        temp_output=$(mktemp)
        sbatch $PROJECT_HOME_DIR/scripts/hpc/run_preprocess.sh "$fastq_file" 2>&1 | tee "$temp_output"
        log $(cat $temp_output)
    else
        log "File not found: $fastq_file"
    fi
done
