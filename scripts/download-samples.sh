#!/bin/bash

PROJECT_HOME_DIR=${PROJECT_HOME_DIR:-/home/projects/22126_NGS/projects/group15}
PROJECT_DATA_DIR=$PROJECT_HOME_DIR/data
LOG_DIR=$PROJECT_HOME_DIR/logs

SRR_IDS=$(cat $PROJECT_HOME_DIR/data/samples.csv | cut -d ',' -f1)

mkdir -p $LOG_DIR
cd $PROJECT_DATA_DIR/raw

for SRR_ID in $SRR_IDS; do
    ftp_url="ftp://$(cat $PROJECT_HOME_DIR/data/filereport_read_run_PRJNA340216.tsv | cut -f7 | grep $SRR_ID)"
    file_name=$(basename $ftp_url)
    if [ -f "$file_name" ]; then
        echo "File already exists, skipping: $file_name"
    else
        echo "Downloading $SRR_ID from $ftp_url ($file_name)"
        wget $ftp_url --append-output $LOG_DIR/SRR-download.log -nv
    fi
done

echo "Done"
