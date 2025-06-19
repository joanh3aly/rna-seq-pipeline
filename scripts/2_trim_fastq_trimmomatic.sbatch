#!/bin/bash
####
#### Joan Healy <joanhealy24@rcsi.ie>
#### Trim FastQ files
#### Inputs: Sample Fastq files from all patients
#### Outputs: Trimmed FastQ files 
####
#SBATCH --job-name=trim_fastq
# #SBATCH --output=/home/data/human_genetics/USERS/jhealy/rna_pipeline_v0.1/scripts/slurmfiles/%j_trim_FASTQ.out
# #SBATCH --error=/home/data/human_genetics/USERS/jhealy/rna_pipeline_v0.1/scripts/slurmfiles/%j_trim_FASTQ.err
#SBATCH --partition=ProdQ
#SBATCH --time=2-00:00:00
#SBATCH --ntasks=1
# # SBATCH --mem=4000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joanhealy24@rcsi.com
# #SBATCH --array=0-1

set -x

# module load BBMap
module load Java/11.0.16 || { echo "ERROR: Failed to load Java/11.0.16 module"; exit 1; }
module load Trimmomatic/0.39-Java-11 || { echo "ERROR: Failed to load Trimmomatic/0.39-Java-11 module"; exit 1; }
# Import parallelize module file
source /home/data/human_genetics/USERS/jhealy/rna_pipeline_v0.1/scripts/parallelize_fastq_module.sbatch



printf "\n\n Trim FASTQ files \n\n"
run_id="$1"
printf "\n\n Run ID TRIM: %s \n\n" "${run_id}"
# WORK_DIR=/home/data/human_genetics/USERS/jhealy/rna_pipeline_v0.1
WORK_DIR="$2"
DELETE_UNPAIRED="${3:-true}" 
INPUT_DIR=$WORK_DIR/$run_id/fastq_inputs/concat_lanes
OUTPUT_DIR_TRIMMED_FASTQ=$WORK_DIR/$run_id/trimmed_fastq
ADAPTER_FASTA="$EBROOTTRIMMOMATIC/adapters/TruSeq3-PE.fa"
TRIMMOMATIC_JAR="$EBROOTTRIMMOMATIC/trimmomatic-0.39.jar"


cleanup_unpaired_files() {
    if [[ "$DELETE_UNPAIRED" == "true" ]]; then
        printf "INFO" "Searching for and deleting unpaired files in '$TRIMMED_FASTQ_DIR'..." >&2
        # Find and delete files whose names contain '_unpaired_'
        find "$TRIMMED_FASTQ_DIR" -maxdepth 1 -type f -name "*_unpaired_*" -delete

        if [ $? -eq 0 ]; then
            printf "INFO" "Unpaired file cleanup completed successfully." >&2
        else
            printf "ERROR" "Failed to clean up some unpaired files. Check permissions or file paths." >&2
        fi
    else
      printf "INFO" "Deletion of unpaired files is disabled by configuration." >&2
    fi
}


# Call the function to populate the fastq_array of input files and their indices from the parallelize_fastqc_module script
if fastq_output=$(populate_fastq_array "$INPUT_DIR"); then
  IFS=$'\n' read -r -d '' -a fastq_arr <<< "$fastq_output"
  unset IFS

  printf "Check FASTQ array in the main script:"
  for index in "${!fastq_arr[@]}"; do
    printf "\n Index: $index, File: ${fastq_arr[$index]}"
  done

  # use the 'fastq_arr' to process fastq files and output fastqc reports
  if [[ -v SLURM_ARRAY_TASK_ID ]]; then
    TASK_ID="$SLURM_ARRAY_TASK_ID"
    printf "\n TASK_ID %s \n" "$TASK_ID"
    if [[ "$TASK_ID" -ge 0 ]] && [[ "$TASK_ID" -lt "${#fastq_arr[@]}" ]]; then
      FASTQ_FILE="${fastq_arr[$TASK_ID]}"
      printf "\n Processing file in array job:%s \n" "$FASTQ_FILE"
      # SAMPLE_NAME=$(basename "$FASTQ_FILE" .fq.gz)
      # printf "\n Sample name: $SAMPLE_NAME" 
      printf "\n Processing: $FASTQ_FILE (Task ID: $TASK_ID)"


      if [[ "$FASTQ_FILE" =~ ^"$INPUT_DIR"/.*1\.fq\.gz$ ]]; then
        printf "\n\n FASTQ_FILE is read 1 (processing only initiates for read 1 & 2 when initial read = 1): %s \n" "${FASTQ_FILE}"
        printf "\n\n FASTQ_FILE read 2 is: %s \n" "${FASTQ_FILE/1.fq.gz/2.fq.gz}"
        # output1="${FASTQ_FILE/.fq.gz/.trim.fq.gz}"
        # printf "\n Output 1: %s \n" "${output1}"
        
        # Extract the filename from the input path
        input_filename=$(basename "$FASTQ_FILE")
        printf "\n Input_filename: %s \n" "${input_filename}"
        # Replace "1.fq.gz" with "1.trim.fq.gz" in the filename
        output_filename_r1="${input_filename/1.fq.gz/1.trim.fq.gz}"
        output_filename_r2="${input_filename/1.fq.gz/2.trim.fq.gz}"
        output_filename_r1_unpaired="${input_filename/1.fq.gz/unpaired_1.trim.fq.gz}"
        output_filename_r2_unpaired="${input_filename/1.fq.gz/unpaired_2.trim.fq.gz}"
        output_filepath_r1="${OUTPUT_DIR_TRIMMED_FASTQ}/${output_filename_r1}"
        output_filepath_r2="${OUTPUT_DIR_TRIMMED_FASTQ}/${output_filename_r2}"
        output_filepath_r1_unpaired="${OUTPUT_DIR_TRIMMED_FASTQ}/${output_filename_r1_unpaired}"
        output_filepath_r2_unpaired="${OUTPUT_DIR_TRIMMED_FASTQ}/${output_filename_r2_unpaired}"
        printf "\n\n Output filepaths: \n R1= %s \n\n R2= %s \n"  "${output_filepath_r1}" "${output_filepath_r2}"

 
        java -jar "$TRIMMOMATIC_JAR" PE \
        -threads 8 \
        -phred33 \
        "${FASTQ_FILE}" \
        "${FASTQ_FILE/1.fq.gz/2.fq.gz}" \
        "${output_filepath_r1}" \
        "${output_filepath_r1_unpaired}" \
      	"${output_filepath_r2}" \
      	"${output_filepath_r2_unpaired}" \
        ILLUMINACLIP:"${ADAPTER_FASTA}":2:30:10 \
        LEADING:3 \
        TRAILING:3 \
        SLIDINGWINDOW:4:15 \
        MINLEN:36 \
        -trimlog "${WORK_DIR}/${run_id}/slurmfiles/${run_id}.trim.log" 

        printf "\n Trimming completed for: $FASTQ_FILE"
        cleanup_unpaired_files
      else
        printf "\n Trimming NOT completed for (read is not 1 / array name doesn't match file name): $FASTQ_FILE"
        # exit 1
      fi
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
printf "\n FastQ Trim processing finished."






