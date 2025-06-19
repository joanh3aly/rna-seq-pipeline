#!/bin/bash
####
#### Joan Healy <joanhealy24@rcsi.ie>
#### FastQC trimmed FastQ files
#### Inputs: Sample trimmed Fastq files from all patients
#### Outputs: FastQC reports
####
#SBATCH --job-name=fastqc_trimmed
# #SBATCH --output=/home/data/human_genetics/USERS/jhealy/rna_pipeline_v0.1/scripts/slurmfiles/%j_fastqc_trimmed.out
# #SBATCH --error=/home/data/human_genetics/USERS/jhealy/rna_pipeline_v0.1/scripts/slurmfiles/%j_fastqc_trimmed.err
#SBATCH --partition=ProdQ
#SBATCH --time=2-00:00:00
#SBATCH --ntasks=1
#SBATCH --mem=4000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joanhealy24@rcsi.com
# # SBATCH --array=0-1%

set -x

module load FastQC
# Parallelize module file
WORK_DIR="$2"
source /home/data/human_genetics/USERS/jhealy/rna_pipeline_v0.1/scripts/parallelize_fastq_module.sbatch
# WORK_DIR=/home/data/human_genetics/USERS/jhealy/rna_pipeline_v0.1

# source ${WORK_DIR}/scripts/parallelize_fastq_module.sbatch

printf "\n\n FASTQC Untrimmed files \n\n"
run_id="$1"
printf "\n\n Run ID : %s \n\n" "${run_id}"
# WORK_DIR="$2"
# # WORK_DIR=/home/data/human_genetics/USERS/jhealy/rna_pipeline_v0.1
# source ${WORK_DIR}/scripts/parallelize_fastq_module.sbatch
INPUT_DIR_FASTQ=$WORK_DIR/$run_id/trimmed_fastq/concat_lanes
OUTPUT_DIR_TRIMMED_FASTQC=$WORK_DIR/$run_id/trimmed_fastqc_outputs


echo "Listing sub-files in: $INPUT_DIR_FASTQ (using find)"
echo "---------------------------------------------------"

# Check if the directory exists
if [[ -d "$INPUT_DIR_FASTQ" ]]; then
  # Find files in the specified directory, with specific options
  # -maxdepth 1: Only look in the specified directory, not its subdirectories
  # -type f: Only find regular files
  # -print: Print the full path of the found items
  find "$INPUT_DIR_FASTQ" -maxdepth 1 -type f -print
else
  echo "Error: Directory '$INPUT_DIR_FASTQ' does not exist." >&2
  exit 1
fi


declare -a fastq_arr=()

# Call the function to populate the fastq_array from the parallelize_fastqc_module script
if fastq_output=$(populate_fastq_array "$INPUT_DIR_FASTQ"); then
#   IFS=$'\n' read -r -d '' -a fastq_arr <<< "$fastq_output"
#   unset IFS
  readarray -t fastq_arr <<< "$fastq_output"

  printf "\n DEBUG: Loop through FASTQ array:"
  for index in "${!fastq_arr[@]}"; do
    printf "\n Index: $index, File: %s \n" "${fastq_arr[$index]}"
  done
  printf "Total files found: %s\n" "${#fastq_arr[@]}"

  # use the 'fastq_arr' to process fastq files and output fastqc reports
  if [[ -v SLURM_ARRAY_TASK_ID ]]; then
    TASK_ID="$SLURM_ARRAY_TASK_ID"
    printf "\n TASK_ID %s \n" "$TASK_ID"
    if [[ "$TASK_ID" -ge 0 ]] && [[ "$TASK_ID" -lt "${#fastq_arr[@]}" ]]; then
      FASTQ_FILE="${fastq_arr[$TASK_ID]}"
      printf "\n Processing file in array job:%s \n" "$FASTQ_FILE"
      SAMPLE_NAME=$(basename "$FASTQ_FILE" .trim.fq.gz)
      printf "\n Sample name: $SAMPLE_NAME" 
      OUTPUT_DIR_DYNAMIC="$OUTPUT_DIR_TRIMMED_FASTQC/${SLURM_JOB_ID}_fastqc_job_${SAMPLE_NAME}" 
      mkdir -p $OUTPUT_DIR_DYNAMIC
      printf "\n Processing: $FASTQ_FILE (Task ID: $TASK_ID)"
       # Run FastQC
    #   fastqc -q "$FASTQ_FILE" --outdir="$OUTPUT_DIR_TRIMMED_FASTQ/${SAMPLE_NAME}_fastqc_job_${SLURM_JOB_ID}"
      fastqc -q "$FASTQ_FILE" --outdir="${OUTPUT_DIR_DYNAMIC}"
      printf "\n FastQC completed for: %s/%s \n" "${OUTPUT_DIR_DYNAMIC}" "$FASTQ_FILE"
      printf "\n FastQC reports available in: %s \n" "$OUTPUT_DIR_TRIMMED_FASTQC"
    else
      printf "\n Error: SLURM_ARRAY_TASK_ID ($TASK_ID) is out of bounds for fastq_arr"
    fi
  else  
   printf "\n Error: SLURM_ARRAY_TASK_ID is not defined (not running as a SLURM array job)"
  fi
else
  printf "\n Failed to populate the FASTQ array." >&2
  exit 1
fi
printf "\n FastQC processing finished."
