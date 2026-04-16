#!/usr/bin/env bash 

# Exit on any error
set -e

# === Input parameters ===

while [[ $# -gt 0 ]]; do
    case "$1" in
        --fastq-dir) 
            fastq_dir="$2"
            shift 2
            ;;
        --output-dir) 
            output_dir="$2"
            shift 2
            ;;
        --fastqc-path) 
            fastqc_path="$2"
            shift 2
            ;;
        --num-threads) 
            num_threads="$2"
            shift 2
            ;;
        --num-parallel) 
            num_parallel="$2"
            shift 2
            ;;
        *) 
            echo "Unknown option: $1"; 
            exit 1;;
    esac
done

# === Validate input parameters ===

if [[ -z "$fastq_dir" || -z "$output_dir" ]]; then
  echo "Usage: $0 --fastq-dir <dir> --output-dir <dir> --fastqc-path <path> [--num-threads <n>] [--num-parallel <n>]"
  exit 1
fi

[[ -z "$num_threads" ]] && num_threads=1
[[ -z "$num_parallel" ]] && num_parallel=1  

# === Read the files in the fastq directory ===

shopt -s nullglob
fastqc_files=("$fastq_dir"/*.fastq.gz "$fastq_dir"/*.fastq)
shopt -u nullglob

# ==== Run the quality check ====

if [[ -z "$fastqc_path" ]]; then
  printf "%s\n" "${fastqc_files[@]}" | parallel -j "$num_parallel" "fastqc -t $num_threads -o $output_dir {}"
else
  echo "Using FastQC from: $fastqc_path"
  printf "%s\n" "${fastqc_files[@]}" | parallel -j "$num_parallel" "$fastqc_path -t $num_threads -o $output_dir {}"
fi
