#!/bin/bash
#set -x
####
#### Joan Healy <joanhealy24@rcsi.ie>
#### Parallelise FastQ files
#### Inputs: Sample FastQ files from all patients
#### Outputs: FastQ file directory array
####
#SBATCH --job-name=parallelize_fastq
#SBATCH --output=/home/data/human_genetics/USERS/jhealy/rna_pipeline_v0.1/scripts/slurmfiles/%j_parallelize_fastq.out
#SBATCH --error=/home/data/human_genetics/USERS/jhealy/rna_pipeline_v0.1/scripts/slurmfiles/%j_parallelize_fastq.err
#SBATCH --partition=ProdQ
#SBATCH --time=2-00:00:00
#SBATCH --ntasks=1
#SBATCH --mem=4000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joanhealy24@rcsi.com
#SBATCH --array=0-1

module load FastQC

# WORK_DIR=/home/data/human_genetics/USERS/jhealy/rna_pipeline_v0.1
# INPUT_DIR=$WORK_DIR/fastq_inputs/urine_batch1/X204SC24042165-Z01-F001/01.RawData/batch1_combined
# WORK_DIR="/home/data/human_genetics/USERS/jhealy/rna_pipeline_v0.1"
# INPUT_DIR="${WORK_DIR}/run_20250522_1017/trimmed_fastq"

echo "Parallelize FastQ module running"
# Function to find and populate an array of FASTQ files
populate_fastq_array() {
  local input_dir="$1"
  declare -a local_fastq_arr=()

  if [[ -d "$input_dir" ]]; then
    for f in "$input_dir"/*; do
      # if [[ "$f" =~ ^"$input_dir"/.*\.fq\.gz$ ]]; then
      if [[ "$f" =~ \.fq\.gz$ ]]; then
        local_fastq_arr+=("$f")
      fi
    done
    # Output the array elements separated by newlines
    printf "%s\n" "${local_fastq_arr[@]}"
  else
    echo "Error: Input directory '$input_dir' does not exist." >&2
    return 1
  fi
}
#   for f in "$input_dir"/*; do
#     # 1. Skip if it's not a regular file (e.g., directories, symlinks, etc.)
#     if [[ ! -f "$f" ]]; then
#         # printf "DEBUG: Skipping non-regular file/directory: %s\n" "$f" >&2 # Uncomment for debugging
#         continue
#     fi

#     # 2. Check if the file ends with .fq.gz or .trim.fq.gz
#     # The regex \.fq\.gz$ will match both sample.fq.gz and sample.trim.fq.gz
#     # So, a single check for the ultimate .fq.gz ending is sufficient if both are desired.
#     if [[ "$f" =~ \.fq\.gz$ ]]; then
#         printf "Found valid FASTQ file for processing: %s\n" "$f" >&2
#         local_fastq_arr+=("$f")
#     else
#         # If it's a file but doesn't end in .fq.gz (e.g., .txt, .log, .R)
#         printf "WARNING: Skipping file with unexpected extension in '%s': %s\n" "$input_dir" "$f" >&2
#     fi
#   done
# }


##############################################################
# TESTING ONLY
# call the function and store the output
# if fastq_output=$(populate_fastq_array "$INPUT_DIR"); then # executes the function named populate_fastq_array and passes the value of the variable $INPUT_DIR as an argument to it. 
#   IFS=$'\n' read -r -d '' -a fastq_array <<< "$fastq_output" # the read command to process the output stored in the fastq_output variable. IFS=$'\n': This temporarily sets the Internal Field Separator (IFS) to a newline character.  IFS is used by Bash to determine word boundaries when splitting strings. # -r: This option prevents backslash escapes from being interpreted. This is important when dealing with file paths that might contain backslashes. -d '':  specifies the delimiter that marks the end of the input. By setting it to an empty string, read will read all the input until the end of the file or the input stream.
#   echo "IFS $IFS"
#   unset IFS

#   echo "Found FASTQ files:"
#   for index in "${!fastq_array[@]}"; do
#     echo "Index: $index, File: ${fastq_array[$index]}" 
#   done
# fi
