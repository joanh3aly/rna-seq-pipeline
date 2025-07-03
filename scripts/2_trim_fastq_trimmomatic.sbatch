#!/bin/bash
####
#### Joan Healy <joanhealy24@rcsi.ie>
#### Trim FastQ files
#### Inputs: Sample Fastq files from all patients
#### Outputs: Trimmed FastQ files 
####
#SBATCH --job-name=trim_fastq
#SBATCH --partition=ProdQ
#SBATCH --time=2-00:00:00
#SBATCH --ntasks=1 --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joanhealy24@rcsi.com
# #SBATCH --array=0-1

set -x -o nounset

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
printf "\n\n DELETE_UNPAIRED?: %s \n\n" "${DELETE_UNPAIRED}"
INPUT_DIR=$WORK_DIR/$run_id/fastq_inputs/concat_lanes
OUTPUT_DIR_TRIMMED_FASTQ=$WORK_DIR/$run_id/trimmed_fastq
ADAPTER_FASTA="$EBROOTTRIMMOMATIC/adapters/TruSeq3-PE.fa"
TRIMMOMATIC_JAR="$EBROOTTRIMMOMATIC/trimmomatic-0.39.jar"


cleanup_unpaired_files_from_scratch () {
    if [[ "$DELETE_UNPAIRED" == "true" ]]; then
        printf "INFO" "Searching for and deleting unpaired files in '$SCRATCH_OUTPUT_DIR'..." >&2
        # Find and delete files whose names contain '_unpaired_'
        find "$SCRATCH_OUTPUT_DIR" -maxdepth 1 -type f -name "*_unpaired_*" -delete

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
      printf "\n Processing: $FASTQ_FILE (Task ID: $TASK_ID)"


      if [[ "$FASTQ_FILE" =~ ^"$INPUT_DIR"/.*1\.fq\.gz$ ]]; then
        printf "\n\n FASTQ_FILE is read 1 (processing only initiates for read 1 & 2 when initial read = 1): %s \n" "${FASTQ_FILE}"
        printf "\n\n FASTQ_FILE read 2 is: %s \n" "${FASTQ_FILE/1.fq.gz/2.fq.gz}"
        
        # Extract the filename from the input path
        input_filename=$(basename "$FASTQ_FILE")
        printf "\n Input_filename: %s \n" "${input_filename}"
        # Replace "1.fq.gz" with "1.trim.fq.gz" in the filename
        input_filepath_NFS_r1="${INPUT_DIR}/${input_filename}"
        input_filepath_NFS_r2="${INPUT_DIR}/${input_filename/1.fq.gz/2.fq.gz}"
        output_filename_r1="${input_filename/1.fq.gz/1.trim.fq.gz}"
        output_filename_r2="${input_filename/1.fq.gz/2.trim.fq.gz}"
        output_filename_r1_unpaired="${input_filename/1.fq.gz/unpaired_1.trim.fq.gz}"
        output_filename_r2_unpaired="${input_filename/1.fq.gz/unpaired_2.trim.fq.gz}"


        # Stage input GDB to fast storage; there is a risk if someone else is also using the scratch space
        # so check for enough space for each file
        #-
        # l_scratch_space=$(df -k --output=avail "${TMPDIR}")
        # l_scratch_space="${l_scratch_space#*Avail*}"
        # echo "Free space on scratch disk (KB):${l_scratch_space}"
        # Get size of the GDB database
        # l_gdb=$(du -ks --apparent-size "${INPUT_DIR}")
        # l_gdb="${l_gdb%[[:blank:]]*}"
        # if (( l_scratch_space > l_gdb )); then
          # File will fit on scratch disk
        echo "Staging input files to scratch disk"
        SCRATCH_INPUT_DIR="${TMPDIR}/trim_input_fastq"
        SCRATCH_OUTPUT_DIR="${TMPDIR}/trim_output_fastq"
        mkdir "${SCRATCH_INPUT_DIR}"
        mkdir "${SCRATCH_OUTPUT_DIR}"

        input_filepath_SCRATCH_r1="$SCRATCH_INPUT_DIR/${input_filename}"
        input_filepath_SCRATCH_r2="$SCRATCH_INPUT_DIR/${input_filename/1.fq.gz/2.fq.gz}"
        cp "${input_filepath_NFS_r1}" "${SCRATCH_INPUT_DIR}"
        cp "${input_filepath_NFS_r2}" "${SCRATCH_INPUT_DIR}"
        
        output_filepath_SCRATCH_r1="${SCRATCH_OUTPUT_DIR}/${output_filename_r1}"
        output_filepath_SCRATCH_r2="${SCRATCH_OUTPUT_DIR}/${output_filename_r2}"
        output_filepath_SCRATCH_r1_unpaired="${SCRATCH_OUTPUT_DIR}/${output_filename_r1_unpaired}"
        output_filepath_SCRATCH_r2_unpaired="${SCRATCH_OUTPUT_DIR}/${output_filename_r2_unpaired}"
 
        java -jar "$TRIMMOMATIC_JAR" PE \
        -threads "$SLURM_CPUS_PER_TASK" \
        -phred33 \
        "${input_filepath_SCRATCH_r1}" \
        "${input_filepath_SCRATCH_r2}" \
        "${output_filepath_SCRATCH_r1}" \
        "${output_filepath_SCRATCH_r1_unpaired}" \
      	"${output_filepath_SCRATCH_r2}" \
      	"${output_filepath_SCRATCH_r2_unpaired}" \
        ILLUMINACLIP:"${ADAPTER_FASTA}":2:30:10 \
        LEADING:3 \
        TRAILING:3 \
        SLIDINGWINDOW:4:15 \
        MINLEN:36 \
        -trimlog "${WORK_DIR}/${run_id}/slurmfiles/${input_filename}.trim.log" 

        printf "\n Trimming completed for: $FASTQ_FILE"
        cleanup_unpaired_files_from_scratch "${SCRATCH_OUTPUT_DIR}"
        cp "${output_filepath_SCRATCH_r1}" "${OUTPUT_DIR_TRIMMED_FASTQ}"
        cp "${output_filepath_SCRATCH_r2}" "${OUTPUT_DIR_TRIMMED_FASTQ}"
        rm -rf "${SCRATCH_INPUT_DIR}" "${SCRATCH_OUTPUT_DIR}"
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
















