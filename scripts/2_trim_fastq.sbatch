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

module load BBMap
# Import parallelize module file
source /home/data/human_genetics/USERS/jhealy/rna_pipeline_v0.1/scripts/parallelize_fastq_module.sbatch

# INPUT_DIR=$WORK_DIR/fastq_inputs/urine_batch1/X204SC24042165-Z01-F001/01.RawData/batch1_combined
# INPUT_DIR=$WORK_DIR/fastq_inputs/urine_batch1/X204SC24042165-Z01-F001/01.RawData/batch1_combined
# OUTPUT_DIR_TRIMMED_FASTQ=$WORK_DIR/trimmed_fastq

printf "\n\n Trim FASTQ files \n\n"
run_id="$1"
printf "\n\n Run ID TRIM: %s \n\n" "${run_id}"
# WORK_DIR=/home/data/human_genetics/USERS/jhealy/rna_pipeline_v0.1
WORK_DIR="$2"
INPUT_DIR=$WORK_DIR/$run_id/fastq_inputs/concat_lanes
OUTPUT_DIR_TRIMMED_FASTQ=$WORK_DIR/$run_id/trimmed_fastq


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
      SAMPLE_NAME=$(basename "$FASTQ_FILE" .fq.gz)
      printf "\n Sample name: $SAMPLE_NAME" 
      printf "\n Processing: $FASTQ_FILE (Task ID: $TASK_ID)"


      if [[ "$FASTQ_FILE" =~ ^"$INPUT_DIR"/.*1\.fq\.gz$ ]]; then
        printf "\n\n FASTQ_FILE is read 1 (processing only initiates for read 1 & 2 when initial read = 1): %s \n" "${FASTQ_FILE}"
        output1="${FASTQ_FILE/.fq.gz/.trim.fq.gz}"
        printf "\n Output 1: %s \n" "${output1}"
        
        # Extract the filename from the input path
        input_filename=$(basename "$FASTQ_FILE")
        printf "\n Input_filename: %s \n" "${input_filename}"
        # Replace "1.fq.gz" with "1.trim.fq.gz" in the filename
        output_filename_r1="${input_filename/1.fq.gz/1.trim.fq.gz}"
        output_filename_r2="${input_filename/1.fq.gz/2.trim.fq.gz}"
        output_filepath_r1="${OUTPUT_DIR_TRIMMED_FASTQ}/${output_filename_r1}"
        output_filepath_r2="${OUTPUT_DIR_TRIMMED_FASTQ}/${output_filename_r2}"
        printf "\n\n Output filepaths: \n R1= %s \n\n R2= %s \n"  "${output_filepath_r1}" "${output_filepath_r2}"

        OUTPUT_DIR_DYNAMIC="$OUTPUT_DIR_TRIMMED_FASTQ/${SLURM_JOB_ID}_stats" 
        # mkdir -p $OUTPUT_DIR_DYNAMIC
        # stats_file="${OUTPUT_DIR_DYNAMIC}/${SLURM_JOB_ID}.stats"
        stats_file="${SLURM_JOB_ID}.stats"
 
        java \
        -Djava.lirary.path=/home/data/tapm/transcriptomics/bbmap/jni/ \
        -ea \
        -Xmx16g \
        -Xms16g \
        -cp /home/data/tapm/transcriptomics/bbmap/current/ \
        jgi.BBDukF \
        -Xmx16g \
        in1=${FASTQ_FILE} \
        in2=${FASTQ_FILE/1.fq.gz/2.fq.gz} \
        out1=${output_filepath_r1} \
        out2=${output_filepath_r2} \
        k=31 \
        mink=5 \
        hdist=1 \
        ktrim=r \
        ref=/home/data/tapm/transcriptomics/bbmap/resources/adapters.fa \
        trimq=20 \
        qtrim=rl \
        stats=$WORK_DIR/$run_id/slurmfiles/${stats_file}.stats \
        threads=8 

        # ./${OUTPUT_DIR_TRIMMED_FASTQ}/${stats_file}.stats \

        printf "\n Trimming completed for: $FASTQ_FILE"
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






