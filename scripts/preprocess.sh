#!/bin/bash
#
# Pre-processing of raw single-end FASTQ reads for alignment
# 1. QC
# 2. Trim adapters, quality-trim and filter
# 3. Post-trim QC
# 4. Decontaminate: Align to human genome and remove unmapped reads

set -e
set -u
set -o pipefail

PROJECT_HOME_DIR=${PROJECT_HOME_DIR:-/home/projects/22126_NGS/projects/group15}
PROJECT_DATA_DIR=$PROJECT_HOME_DIR/data
PROEJCT_SCRIPTS_DIR=$PROJECT_HOME_DIR/scripts
TRIMMED_DIR=$PROJECT_DATA_DIR/trimmed
DECONTAMINATED_DIR=$PROJECT_DATA_DIR/decontaminated
LOG_DIR=$PROJECT_HOME_DIR/logs

threads=4

mkdir -p -m 775 $LOG_DIR
mkdir -p -m 775 $TRIMMED_DIR
mkdir -p -m 775 $DECONTAMINATED_DIR

log_file="$LOG_DIR/$(basename $0).log"


log() {
    local msg="[$(date +"%Y-%m-%d %H:%M:%S")] [$(basename $0)] $*" 
    echo "$msg" >> "$log_file" 2>&1
    echo "$msg"
}

qc() {
    $PROEJCT_SCRIPTS_DIR/qc-samples.sh "$1"
}

trim() {
    local input_fastq=$1
    local trimmed_fastq=$2
    local trimmomatic="java -jar /home/ctools/Trimmomatic-0.39/trimmomatic-0.39.jar"
    local adapters="/home/ctools/Trimmomatic-0.39/adapters/NexteraPE-PE.fa"

    if [ -f "$trimmed_fastq" ]; then
        log "Trimmed file already exists ($trimmed_fastq), skipping $input_fastq"
    else
        log "Trimming $input_fastq, output: $trimmed_fastq"
        temp_output=$(mktemp)

        $trimmomatic SE \
        -threads $threads \
        -phred33 \
        "$input_fastq" \
        "$trimmed_fastq" \
        ILLUMINACLIP:$adapters:2:30:10 \
        LEADING:20 \
        TRAILING:20 \
        SLIDINGWINDOW:4:20 \
        MINLEN:60 \
        2>&1 | tee "$temp_output"

        log $(cat $temp_output)
        log "Done trimming $input_fastq, output: $trimmed_fastq"
    fi
}

remove_human_contaminants() {
    local input_fastq=$1
    local output_fastq=$2
    local human_ref=/home/databases/references/human/GRCh38_full_analysis_set_plus_decoy_hla.fa

    if [ -f "$output_fastq" ]; then
        log "Decontaminated file already exists ($output_fastq), skipping $input_fastq"
    else
        log "Decontaminating $input_fastq reads: aligning to $human_ref and removing unmapped reads"
        bwa mem -t "$threads" "$human_ref" "$input_fastq" \
            | samtools view -b -f 4 \
            | samtools fastq \
            | gzip > "$output_fastq"
        log "Done decontaminating $input_fastq, output: $output_fastq"
    fi
}


if [ -z "$1" ]; then
    echo "Usage:
    bash preprocess.sh myfile.fastq.gz
or
    bash preprocess.sh mydir/"
    exit 1
fi

if [ -f $1 ]; then
    input=$1
elif [ -d $1 ]; then
    input=$1/*.fastq.gz
else
    log "Invalid input: must either be a directory or a single file"
    exit 1
fi

shopt -s nullglob
for fastq_file in $input; do
    log "==== Pre-processing $fastq_file ===="
    sample_name=$(echo $(basename $fastq_file) | cut -d '.' -f1)
    trimmed_fastq="$TRIMMED_DIR/${sample_name}_trimmed.fastq.gz"
    decontaminated_fastq="$DECONTAMINATED_DIR/${sample_name}_trimmed_nonhuman.fastq.gz"

    qc "$fastq_file"
    trim "$fastq_file" "$trimmed_fastq"
    qc "$trimmed_fastq"
    remove_human_contaminants "$trimmed_fastq" "$decontaminated_fastq"
    log "==== Done pre-processing $fastq_file ===="
done
