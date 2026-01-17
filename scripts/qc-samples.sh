#!/bin/bash

set -e
set -u
set -o pipefail

PROJECT_HOME_DIR=${PROJECT_HOME_DIR:-/home/projects/22126_NGS/projects/group15}
PROJECT_DATA_DIR=$PROJECT_HOME_DIR/data
LOG_DIR=$PROJECT_HOME_DIR/logs
PROJECT_DATA_FASTQ_DIR=$PROJECT_DATA_DIR/fastqc

log_file="$LOG_DIR/$(basename $0).log"
fastqc=/home/ctools/FastQC/fastqc

log() {
    local msg="[$(date +"%Y-%m-%d %H:%M:%S")] [$(basename $0)] $*" 
    echo "$msg" >> "$log_file" 2>&1
    echo "$msg"
}

mkdir -p $LOG_DIR


if [ -z "$1" ]; then
    echo "Usage:
    bash qc-samples.sh myfile.fastq.gz
or
    bash qc-samples.sh mydir/"
    exit 1
fi

if [ -f $1 ]; then
    input=$1
elif [ -d $1 ]; then
    input=$1/*.fastq.gz
else
    log "Invalid input $1: must either be a directory or a single file"
    exit 1
fi

shopt -s nullglob
for fastq_file in $input; do
    file_name="$(basename "$fastq_file")"
    qc_report="$PROJECT_DATA_FASTQ_DIR"/"$(echo "$file_name" | cut -d '.' -f1)_fastqc.html"
    if [ -f "$qc_report" ]; then
        log "QC report already exists ($qc_report), skipping $fastq_file"
    else
        log "Running FastQC on $fastq_file"
        temp_output=$(mktemp)
        $fastqc -o "$PROJECT_DATA_FASTQ_DIR" "$fastq_file" 2>&1 | tee "$temp_output"
        log $(cat $temp_output)
        log "Done running FastQC on $fastq_file"
    fi
done
