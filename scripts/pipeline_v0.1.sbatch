#!/bin/bash
####
#### Joan Healy <joanhealy24@rcsi.ie>
#### Main Pipeline V0.1
#### Inputs: Base directory
#### Outputs: Entire pipeline
####
#SBATCH --job-name=pipeline_v0.1
#SBATCH --output=/home/data/human_genetics/USERS/jhealy/rna_pipeline_v0.1/scripts/slurmfiles/%j_pipeline_v0.1.out
#SBATCH --error=/home/data/human_genetics/USERS/jhealy/rna_pipeline_v0.1/scripts/slurmfiles/%j_pipeline_v0.1.err
#SBATCH --partition=ProdQ
#SBATCH --time=2-00:00:00
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joanhealy24@rcsi.com

set -x 
module load R
WORK_DIR=/home/data/human_genetics/USERS/jhealy/rna_pipeline_v0.1
INPUT_DROP=$WORK_DIR/fastq_drop
# INDEX="${WORK_DIR}/reference_data/human/GRCh38/HISAT2_index/p14_masked/GRCh38_basic_ss_exon_ind"
INDEX="${WORK_DIR}/reference_data/human/GRCh38/HISAT2_index/p14_masked_gtf_masked/GRCh38_basic_ss_exon_index"
# INDEX="${WORK_DIR}/reference_data/human/GRCh38/HISAT2_index/prebuilt_grch38_snp_tran/genome_snp_tran"

# R_LIBRARIES_PATH="${OUTPUT_DIR}/libraries/"


generate_run_id_directory() {
    local current_date=$(date +%Y%m%d)
    local current_time=$(date +%H%M)
    # Construct the run identifier
    local run_id="run_${current_date}_${current_time}"
    local output_directory="${WORK_DIR}/${run_id}" # Removed trailing slash here for consistency, add it back if needed later

    # Create the main run directory
    mkdir -p "${output_directory}" || { echo "ERROR: Failed to create main run directory: ${output_directory}" >&2; return 1; }
    printf "Created output directory: %s\n" "${output_directory}\n" >&2 # Send to stderr, added newline for clean debug output

    # Create subdirectories
    mkdir -p "${output_directory}/fastq_inputs" || { echo "ERROR: Failed to create fastq_inputs directory." >&2; return 1; }
    mkdir -p "${output_directory}/fastq_inputs/concat_lanes" || { echo "ERROR: Failed to create concat_lanes directory." >&2; return 1; }
    mkdir -p "${output_directory}/aligner_outputs" || { echo "ERROR: Failed to create aligner_outputs directory." >&2; return 1; }
    mkdir -p "${output_directory}/slurmfiles" || { echo "ERROR: Failed to create slurmfiles directory." >&2; return 1; }
    mkdir -p "${output_directory}/r_outputs" || { echo "ERROR: Failed to create r_outputs directory." >&2; return 1; }
    mkdir -p "${output_directory}/trimmed_fastq" || { echo "ERROR: Failed to create trimmed_fastq directory." >&2; return 1; }
    mkdir -p "${output_directory}/trimmed_fastqc_outputs" || { echo "ERROR: Failed to create trimmed_fastqc_outputs directory." >&2; return 1; }
    mkdir -p "${output_directory}/untrimmed_fastqc_outputs" || { echo "ERROR: Failed to create untrimmed_fastqc_outputs directory." >&2; return 1; }
    mkdir -p "${output_directory}/psiclass_outputs" || { echo "ERROR: Failed to create psiclass_outputs directory." >&2; return 1; }

    echo "DEBUG (in generate_run_id_directory): INPUT_DROP is: ${INPUT_DROP}" >&2

    # Copy all files ending in .fq.gz under the fastq_drop directory
    # Note: Assumed INPUT_DROP is WORK_DIR/fastq_drop
    find "${INPUT_DROP}" -type f -name "*.fq.gz" -print0 | while IFS= read -r -d $'\0' file_to_copy; do
        local base_name=$(basename "${file_to_copy}")
        echo "DEBUG: Trying to copy file: '${base_name}' to '${output_directory}/fastq_inputs'" >&2
        cp "${file_to_copy}" "${output_directory}/fastq_inputs/${base_name}"
        if [ $? -ne 0 ]; then
            echo "ERROR: Failed to copy '${file_to_copy}' to '${output_directory}/fastq_inputs/'" >&2
            # Decide if you want to exit the entire script here or just log the error and continue
            # For now, it just logs, but a 'return 1' here would stop the function on first copy error.
        fi
    done

    # This is the ONLY line that should output to stdout
    echo "${output_directory}" # No trailing slash for easier variable use, add it if needed for specific commands.
}


concatenate_sample_lanes() {
    local RAW_FASTQ_DIR="$1/batch1+2_combined"
    local COMBINED_FASTQ_DIR="$3/$2/fastq_inputs/concat_lanes"

    echo "Starting FASTQ lane combination process..."
    echo "Raw FASTQ directory: $RAW_FASTQ_DIR"
    echo "Combined FASTQ directory: $COMBINED_FASTQ_DIR"
    echo ""

  #  readarray -t unique_samples < <(find "$RAW_FASTQ_DIR" -maxdepth 1 -name "*.fq.gz" -o -name "*.fastq.gz" | \
  #                                    sed -E 's|.*/||; s/_L[0-9]+_[12]\.(fq|fastq)\.gz$//' | \
  #                                    sort -u)

  
    readarray -t unique_samples < <(find "$RAW_FASTQ_DIR" -maxdepth 1 -name "*.fq.gz" -o -name "*.fastq.gz" | \
                                     sed -E 's|.*/||' | \
                                     awk -F_ '{print $1}' | \
                                     sort -u)

    # --- Check if any samples were found ---
    if [ ${#unique_samples[@]} -eq 0 ]; then
        echo "ERROR: No FASTQ files found in '$RAW_FASTQ_DIR' or pattern doesn't match the expected naming convention ('*_L<lane>_<read_pair>.<fq/fastq>.gz')."
        return 1 # Return a non-zero status to indicate an error
    fi

    echo "Found unique samples:"
    for sample_id in "${unique_samples[@]}"; do
        echo "- $sample_id"
    done
    echo ""

    # --- Step 2: Iterate through each unique sample and combine its lanes ---
    for sample_id in "${unique_samples[@]}"; do
        echo "Processing sample: $sample_id"

        # --- Combine Read 1 files ---
        # Need to specifically match .fq.gz or .fastq.gz
        input_r1_files_patterns=("$RAW_FASTQ_DIR/${sample_id}_L*_1.fq.gz" "$RAW_FASTQ_DIR/${sample_id}_L*_1.fastq.gz")
        output_r1_file="${COMBINED_FASTQ_DIR}/${sample_id}_1.fq.gz" # Outputting gzipped FASTQ

        declare -a r1_files=()
        # Find all matching files (both .fq.gz and .fastq.gz)
        while IFS= read -r -d $'\0' file; do
            r1_files+=("$file")
        done < <(find "$RAW_FASTQ_DIR" -maxdepth 1 \
                    -name "${sample_id}_*_1.fq.gz" -print0 -o \
                    -name "${sample_id}_*_1.fastq.gz" -print0)

        # done < <(find "$RAW_FASTQ_DIR" -maxdepth 1 -name "${sample_id}_L*_1.fq.gz" -print0 -o -name "${sample_id}_L*_1.fastq.gz" -print0)


        if [ ${#r1_files[@]} -gt 0 ]; then
            if [ -f "$output_r1_file" ]; then
                echo "  WARNING: Output R1 file '$output_r1_file' already exists. Skipping combining for R1."
            else
                echo "  Combining R1 files for '$sample_id' (input is gzipped)..."
                # Use zcat to decompress input files before piping to gzip for re-compression
                zcat "${r1_files[@]}" | gzip > "$output_r1_file"
                echo "  Combined R1 into: $output_r1_file"
            fi
        else
            echo "  No R1 files found for '$sample_id'. Skipping R1 combining."
        fi

        # --- Combine Read 2 files (for paired-end) ---
        input_r2_files_patterns=("$RAW_FASTQ_DIR/${sample_id}_L*_2.fq.gz" "$RAW_FASTQ_DIR/${sample_id}_L*_2.fastq.gz")
        output_r2_file="${COMBINED_FASTQ_DIR}/${sample_id}_2.fq.gz"

        declare -a r2_files=()
        while IFS= read -r -d $'\0' file; do
            r2_files+=("$file")
        done < <(find "$RAW_FASTQ_DIR" -maxdepth 1 \
                    -name "${sample_id}_*_2.fq.gz" -print0 -o \
                    -name "${sample_id}_*_2.fastq.gz" -print0) 

        # done < <(find "$RAW_FASTQ_DIR" -maxdepth 1 -name "${sample_id}_L*_2.fq.gz" -print0 -o -name "${sample_id}_L*_2.fastq.gz" -print0)

        if [ ${#r2_files[@]} -gt 0 ]; then
            if [ -f "$output_r2_file" ]; then
                echo "  WARNING: Output R2 file '$output_r2_file' already exists. Skipping combining for R2."
            else
                echo "  Combining R2 files for '$sample_id' (input is gzipped)..."
                zcat "${r2_files[@]}" | gzip > "$output_r2_file"
                echo "  Combined R2 into: $output_r2_file"
            fi
        else
            echo "  No R2 files found for '$sample_id'. Assuming single-end data if R1 was found."
        fi

        echo ""
    done

    echo "FASTQ lane combination completed for all samples."
    return 0 # Return success status
}



# Generate RUN ID, directory structure, and move all input files from fastq_drop to the run ID dir's fastq_inputs dir ######################################
echo "Starting dir generation..."
STAGE1_START=$SECONDS

echo "DEBUG: Before calling generate_run_id_directory()"
# Generate the RUN_ID and output directory once
run_output_dir="$(generate_run_id_directory)"
echo "DEBUG: After calling generate_run_id_directory()"
echo "DEBUG: RUN_OUTPUT_DIR is: $run_output_dir"
RUN_ID=$(basename "$run_output_dir") # Extract just the run_id part

STAGE1_DURATION=$(( SECONDS - STAGE1_START ))
echo "Finished dir generation in $STAGE1_DURATION seconds."


# --- Main execution block ---
# if [ "$#" -ne 2 ]; then
#     echo "Usage: $0 <raw_fastq_directory> <combined_fastq_directory>"
#     echo "  Example: $0 ./data/raw_fastqs ./data/combined_fastqs"
#     exit 1
# fi


# Concatenate samples ################################################################
echo "Starting sample concatenation..."
STAGE2_START=$SECONDS

concatenate_sample_lanes "$INPUT_DROP" "$RUN_ID" "$WORK_DIR"

STAGE2_DURATION=$(( SECONDS - STAGE2_START ))
echo "Finished sample concatenation in $STAGE2_DURATION seconds."

# Parallelisation job array arguemnts ###########################################################################
echo "Starting parallelisation..."
STAGE3_START=$SECONDS
# Define the input FASTQ directory for the input files for the *current* run
CURRENT_RUN_FASTQ_DIR="${WORK_DIR}/${RUN_ID}/fastq_inputs/concat_lanes"
echo "DEBUG CURRENT_RUN_FASTQ_DIR is: ${CURRENT_RUN_FASTQ_DIR}"

# Get the actual number of FASTQ files copied into the run-specific directory
NUM_FASTQ_FILES=$(find "$CURRENT_RUN_FASTQ_DIR" -maxdepth 1 -type f -name "*.fq.gz" | wc -l)
printf "\n\nDEBUG NUM_FASTQ_FILES is: ${NUM_FASTQ_FILES} \n\n"
#  0-indexed array tasks
MAX_ARRAY_INDEX=$((NUM_FASTQ_FILES - 1))
printf "\n\nDEBUG MAX_ARRAY_INDEX is: ${MAX_ARRAY_INDEX} \n\n"

STAGE3_DURATION=$(( SECONDS - STAGE3_START ))
echo "Finished parallelisation in $STAGE3_DURATION seconds."

# FastQC untrimmed files ###########################################################################
echo "Starting FASTQC untrimmed..."
STAGE4_START=$SECONDS

echo "Submitting FastQC array job for $NUM_FASTQ_FILES files (0-${MAX_ARRAY_INDEX})" #with concurrency $MAX_CONCURRENT_TASKS."
sbatch --job-name=fastqc_untrimmed \
       --array="0-${MAX_ARRAY_INDEX}" \
       --output="/${WORK_DIR}/${RUN_ID}/slurmfiles/%A_%a_fastqc_untrimmed.out" \
       --error="/${WORK_DIR}/${RUN_ID}/slurmfiles/%A_%a_fastqc_untrimmed.err" \
       ./1_fastqc_untrimmed.sbatch "$RUN_ID" "$WORK_DIR"

STAGE4_DURATION=$(( SECONDS - STAGE4_START ))
echo "Finished FastQC untrimmed in $STAGE4_DURATION seconds."

# FASTQ trim files with BBMAP ######################################################################################
# echo "Submitting trim FastQ array job for $NUM_FASTQ_FILES files (0-${MAX_ARRAY_INDEX})" #with concurrency $MAX_CONCURRENT_TASKS."
# TRIMMING_JOB_ID=$(sbatch --parsable \
#        --job-name=trim_fastq \
#        --mem=20G \
#        --array="0-${MAX_ARRAY_INDEX}" \
#        --output="/${WORK_DIR}/${RUN_ID}/slurmfiles/%A_%a_trim_FASTQ.out" \
#        --error="/${WORK_DIR}/${RUN_ID}/slurmfiles/%A_%a_trim_FASTQ.err" \
#        ./2_trim_fastq_bbmap.sbatch "$RUN_ID" "$WORK_DIR")

# if [ -z "$TRIMMING_JOB_ID" ]; then
#     echo "ERROR: Trimming job submission failed." >&2
#     exit 1
# fi
# echo "Trimming job ID: $TRIMMING_JOB_ID"

# FASTQ trim files with TRIMMOMATIC ######################################################################################
echo "Starting FASTQ trimming..."
STAGE5_START=$SECONDS

echo "Submitting trim FastQ array job for $NUM_FASTQ_FILES files (0-${MAX_ARRAY_INDEX})" #with concurrency $MAX_CONCURRENT_TASKS."
TRIMMING_JOB_ID=$(sbatch --parsable \
       --job-name=trim_fastq \
       --mem=200G \
       --array="0-${MAX_ARRAY_INDEX}" \
       --output="/${WORK_DIR}/${RUN_ID}/slurmfiles/%A_%a_trim_FASTQ.out" \
       --error="/${WORK_DIR}/${RUN_ID}/slurmfiles/%A_%a_trim_FASTQ.err" \
       ./2_trim_fastq_trimmomatic.sbatch "$RUN_ID" "$WORK_DIR" "true")

if [ -z "$TRIMMING_JOB_ID" ]; then
    echo "ERROR: Trimming job submission failed." >&2
    exit 1
fi
echo "Trimming job ID: $TRIMMING_JOB_ID"

STAGE5_DURATION=$(( SECONDS - STAGE5_START ))
echo "Finished trimming in $STAGE5_DURATION seconds."

# FASTQC trimmed files #####################################################################################
echo "Starting FASTQC of trimmed FastQ..."
STAGE6_START=$SECONDS

echo "Submitting trimmed FastQC array job for $NUM_FASTQ_FILES files (0-${MAX_ARRAY_INDEX})" #with concurrency $MAX_CONCURRENT_TASKS."
FASTQC_TRIMMED_JOB_ID=$(sbatch --parsable \
       --job-name=fastqc_trimmed \
       --array="0-${MAX_ARRAY_INDEX}" \
       --output="/${WORK_DIR}/${RUN_ID}/slurmfiles/%A_%a_fastqc_trimmed.out" \
       --error="/${WORK_DIR}/${RUN_ID}/slurmfiles/%A_%a_fastqc_trimmed.err" \
       --dependency=afterany:$TRIMMING_JOB_ID \
       ./3_fastqc_trimmed.sbatch "$RUN_ID" "$WORK_DIR")

if [ -z "$FASTQC_TRIMMED_JOB_ID" ]; then
    echo "ERROR: Trimmed FASTQC job submission failed." >&2
    exit 1
fi
echo "Trimmed FASTQC job ID: $FASTQC_TRIMMED_JOB_ID"

STAGE6_DURATION=$(( SECONDS - STAGE6_START ))
echo "Finished trimmed FastQC in $STAGE6_DURATION seconds."

# ALIGN FASTQ trimmed files #####################################################################################
echo "Starting alignment..."

printf "\nWaiting for trimming job (%s) to complete...\n" "$TRIMMING_JOB_ID"
# Loop until the job is no longer in the queue
while squeue -h -j "$TRIMMING_JOB_ID" &>/dev/null; do
    echo "Trimming job $TRIMMING_JOB_ID is still running. Waiting 30 seconds..."
    sleep 30
done
# Check final status using sacct to ensure it completed successfully
JOB_STATE=$(sacct -j "$TRIMMING_JOB_ID" --format=State -n -P | head -n1)
JOB_STATE=$(echo "$JOB_STATE" | tr -d '[:space:]') 

if [[ "$JOB_STATE" == "COMPLETED" ]]; then
    echo "Trimming job $TRIMMING_JOB_ID finished successfully. Proceeding to calculate alignment index."
    TRIMMED_FASTQ_DIR="${WORK_DIR}/${RUN_ID}/trimmed_fastq"
    # Get the actual number of trimmed FASTQ files copied into the trimmed fastq directory
    NUM_TRIMMED_FASTQ_FILES=$(find "$TRIMMED_FASTQ_DIR" -maxdepth 1 -type f -name "*.fq.gz" | wc -l)
    printf "\n\n DEBUG: NUM_TRIMMED_FASTQ_FILES is: ${NUM_TRIMMED_FASTQ_FILES} \n\n"
    #  0-indexed array tasks
    MAX_TRIMMED_ARRAY_INDEX=$((NUM_TRIMMED_FASTQ_FILES - 1))
    printf "\n\n DEBUG: MAX_TRIMMED_ARRAY_INDEX is: ${MAX_TRIMMED_ARRAY_INDEX} \n\n"

    echo "Submitting alignment FastQ array job for $NUM_TRIMMED_FASTQ_FILES files (0-${MAX_TRIMMED_ARRAY_INDEX})" 
    sbatch --parsable \
       --job-name=hisat2_align \
       --array="0-${MAX_TRIMMED_ARRAY_INDEX}" \
       --output="${WORK_DIR}/${RUN_ID}/slurmfiles/%A_%a_hisat2_align_trimmed.out" \
       --error="${WORK_DIR}/${RUN_ID}/slurmfiles/%A_%a_hisat2_align_trimmed.err" \
       --dependency=afterany:$TRIMMING_JOB_ID \
       ./4_hisat2_align_multiple.sbatch "$RUN_ID" "$WORK_DIR" "$INDEX"
else
    echo "ERROR: Trimming job $TRIMMING_JOB_ID finished with state: $JOB_STATE. Alignment will not proceed." >&2
    exit 1 # Exit if the trimming job failed
fi


echo "Submitting alignment FastQ array job for $NUM_TRIMMED_FASTQ_FILES files (0-${MAX_TRIMMED_ARRAY_INDEX})" #with concurrency $MAX_CONCURRENT_TASKS."
# HISAT2_ALIGN_TRIMMED_JOB_ID=$(sbatch --parsable \
#        --job-name=hisat2_align \
#        --array="0-${MAX_TRIMMED_ARRAY_INDEX}" \
#        --output="/${WORK_DIR}/${RUN_ID}/slurmfiles/%A_%a_hisat2_align_trimmed.out" \
#        --error="/${WORK_DIR}/${RUN_ID}/slurmfiles/%A_%a_hisat2_align_trimmed.err" \
#        --dependency=afterany:$TRIMMING_JOB_ID \
#        ./4_hisat2_align_multiple.sbatch "$RUN_ID" "$WORK_DIR" "$INDEX")

# if [ -z "$HISAT2_ALIGN_TRIMMED_JOB_ID" ]; then
#     echo "ERROR: Alignment job submission failed." >&2
#     exit 1
# fi
# echo "HISAT2 Alignment job ID: $HISAT2_ALIGN_TRIMMED_JOB_ID"



#######################################################################################
# Submit R script, passing the variable as a command-line argument
# WORK_DIR=/home/data/human_genetics/USERS/jhealy/rna_pipeline_v0.1
# RUN_ID=run_20250522_1017
# # R_LIBRARY_PATH=
# Rscript 4_index.R  "$WORK_DIR" "$RUN_ID" "$R_LIBRARY_PATH"

# if ${index-genome} is true:
#   sbatch 4_index_ref.sbatch "${index-genome}" "${aligner-type}"

#