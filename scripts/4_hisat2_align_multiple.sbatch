#!/bin/bash
####
#### Joan Healy <joanhealy24@rcsi.ie>
#### Main Pipeline V0.1
#### Inputs: Base directory
#### Outputs: Entire pipeline
####
#SBATCH --job-name=hisat2_align
#SBATCH --output=hisat2_align_%A_%a.out
#SBATCH --error=hisat2_align_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16     
#SBATCH --mem=100G              
#SBATCH --time=08:00:00       
#SBATCH --partition=ProdQ 
#SBATCH --array=0-1

set -x

# Load modules 
module load Python/3.11.5-GCCcore-13.2.0
module load HISAT2/2.2.1-foss-2019b
module load SAMtools/1.19.2-GCC-13.2.0
WORK_DIR="$2"
# WORK_DIR=/home/data/human_genetics/USERS/jhealy/rna_pipeline_v0.1
source $WORK_DIR/scripts/parallelize_fastq_module.sbatch

RUN_ID="$1"
# RUN_ID="run_20250620_2125"
# DELETE_UNPAIRED="${3:-true}" 

TRIMMED_FASTQ_DIR="${WORK_DIR}/$RUN_ID/trimmed_fastq"
SLURM_SUBMIT_DIR="${WORK_DIR}/scripts"

HISAT2_INDEX="$3"
# HISAT2_INDEX="${WORK_DIR}/reference_data/human/GRCh38/HISAT2_index/v_2.2.1-foss-2019b/GRCh38_splices_index"
# HISAT2_INDEX="${WORK_DIR}/reference_data/human/GRCh38/grch38_snp/genome_snp
# HISAT2_INDEX="${WORK_DIR}/reference_data/human/GRCh38/HISAT2_index/prebuilt_grch38_snp_tran/genome_snp_tran"
# HISAT2_INDEX="${WORK_DIR}/reference_data/human/GRCh38/HISAT2_index/v_2.2.1-foss-2019b_basic"
# HISAT2_INDEX="${WORK_DIR}/reference_data/human/GRCh38/HISAT2_index/v_2.2.1-foss-2019b_masked/GRCh38_splices_index_masked"
# HISAT2_INDEX="${WORK_DIR}/reference_data/human/GRCh38/HISAT2_index/v_2.2.1-foss-2019b_pc_trans_masked/GRCh38_pc_trans_splices_index_masked"

OUTPUT_DIR="${WORK_DIR}/$RUN_ID/aligner_outputs"


# Call the function to populate the fastq_array of input files and their indices from the parallelize_fastqc_module script
if fastq_output=$(populate_fastq_array "$TRIMMED_FASTQ_DIR"); then
  IFS=$'\n' read -r -d '' -a fastq_arr <<< "$fastq_output"
  unset IFS

  printf "DEBUG: Check FASTQ array in the main script..."
  for index in "${!fastq_arr[@]}"; do
    printf "\n Index: $index, File: ${fastq_arr[$index]}"
  done

  # use the 'fastq_arr' to process fastq files and output fastqc reports
  if [[ -v SLURM_ARRAY_TASK_ID ]]; then
    TASK_ID="$SLURM_ARRAY_TASK_ID"
    printf "\n TASK_ID %s \n" "$TASK_ID"
    if [[ "$TASK_ID" -ge 0 ]] && [[ "$TASK_ID" -lt "${#fastq_arr[@]}" ]]; then
      FASTQ_ARRAY_FILE="${fastq_arr[$TASK_ID]}"
      printf "\n Processing file in array job:%s \n" "$FASTQ_ARRAY_FILE"
      SAMPLE_NAME=$(basename "$FASTQ_ARRAY_FILE" _1.trim.fq.gz)
      printf "\n Sample name: $SAMPLE_NAME" 
      printf "\n Processing: $FASTQ_FILE (Task ID: $TASK_ID)"


      # if [[ "$FASTQ_ARRAY_FILE" =~ ^"$TRIMMED_FASTQ_DIR"/.*1\.trim\.fq\.gz$ ]]; then
      if [[ "$FASTQ_ARRAY_FILE" =~ ^"$TRIMMED_FASTQ_DIR"/.*_1\.trim\.fq\.gz$ ]] && \
        [[ ! "$FASTQ_ARRAY_FILE" =~ _unpaired_ ]]; then
        printf "\n\n FASTQ_FILE is a valid paired Read 1: %s \n" "${FASTQ_ARRAY_FILE}"
        
        # Extract the filename from the input path
        INPUT_FILENAME_R1=$(basename "$FASTQ_ARRAY_FILE")
        INPUT_FILENAME_R2="${INPUT_FILENAME_R1/_1.trim.fq.gz/_2.trim.fq.gz}"
        printf "\n Input_filename R1: %s \n\n Input_filename R2: %s \n" "${INPUT_FILENAME_R1}" "${INPUT_FILENAME_R2}"
        INPUT_FILEPATH_R1="${TRIMMED_FASTQ_DIR}/${INPUT_FILENAME_R1}"
        INPUT_FILEPATH_R2="${TRIMMED_FASTQ_DIR}/${INPUT_FILENAME_R2}"

        if [[ ! -f "$INPUT_FILEPATH_R2" ]]; then
          printf "\n WARNING: R2 file (%s) corresponding to R1 (%s) not found. Skipping this pair.\n" \
              "${INPUT_FILEPATH_R2}" "${INPUT_FILEPATH_R1}"
          continue 
        fi
        printf "\n\n INput filepaths: \n R1= %s \n\n R2= %s \n"  "${INPUT_FILEPATH_R1}" "${INPUT_FILEPATH_R2}"
        OUTPUT_SAM="${OUTPUT_DIR}/${SAMPLE_NAME}.sam"
        OUTPUT_BAM="${OUTPUT_DIR}/${SAMPLE_NAME}.bam"
        OUTPUT_SORTED_BAM="${OUTPUT_DIR}/${SAMPLE_NAME}.sorted.bam"
        OUTPUT_INDEX_BAI="${OUTPUT_DIR}/${SAMPLE_NAME}.sorted.bam.bai"
        ALIGN_LOG="${OUTPUT_DIR}/${SAMPLE_NAME}.hisat2.log"

        echo "Starting HISAT2 alignment for ${INPUT_FILENAME_R1}..."
        echo "HISAT2 Index: ${HISAT2_INDEX}"
        echo "Output Sorted BAM: ${OUTPUT_SORTED_BAM}"

        hisat2 -x "$HISAT2_INDEX" \
                -1 "$INPUT_FILEPATH_R1" \
                -2 "$INPUT_FILEPATH_R2" \
                -p "$SLURM_CPUS_PER_TASK" \
                --summary-file "$ALIGN_LOG" \
                --avoid-pseudogene \
                --dta \
                | samtools view -bS - > "$OUTPUT_BAM"

        printf "\n Alignment completed for: $FASTQ_ARRAY_FILE"
      else
        printf "\n Alignment NOT completed for read 2 file - reverse reads are handled by paired read 1 processing"
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


# Check if alignment was successful, then create BAM indexing file.
if [ $? -eq 0 ]; then
    echo "\n\n HISAT2 alignment complete for ${SAMPLE_NAME}. Sorting BAM file..."
    samtools sort -@ "${SLURM_CPUS_PER_TASK}" \
              "${OUTPUT_BAM}" \
              -o "${OUTPUT_SORTED_BAM}"
    if [ $? -eq 0 ]; then
        echo "BAM file sorted. Indexing..."
        samtools index "${OUTPUT_SORTED_BAM}" \
                "${OUTPUT_INDEX_BAI}"
        if [ $? -eq 0 ]; then
            # echo "BAM file indexed: ${OUTPUT_INDEX_BAI}"
            echo "Alignment process completed successfully for ${SAMPLE_NAME}."
        else
            echo "Error: Failed to index BAM file for ${SAMPLE_NAME}."
            exit 1
        fi
    else
        echo "Error: Failed to sort BAM file for ${SAMPLE_NAME}."
        exit 1
    fi
else
    echo "Error: HISAT2 alignment failed for ${SAMPLE_NAME}."
    exit 1
fi

# Optional: Remove unsorted BAM to save space (uncomment if desired)
# rm "$OUTPUT_BAM"

echo "Script finished for ${SAMPLE_NAME}."