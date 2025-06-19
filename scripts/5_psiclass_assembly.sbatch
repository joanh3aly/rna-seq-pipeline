#!/bin/bash
#
# SLURM script to run PsiCLASS transcriptome assembly on HISAT2 BAM files.
#
# Usage:
#   sbatch 5_psiclass_assembly.sbatch <reference_fasta> <reference_gtf> <output_dir> <bam_list_file>
#
# Example (assuming you have a file named 'bam_files.txt' with paths to your BAMs):
#   sbatch 5_psiclass_assembly.sbatch \
#     /path/to/GRCh38.p14.genome.fa \
#     /path/to/gencode.v44.annotation.gtf \
#     /path/to/psiclass_output_dir \
#     /path/to/bam_files.txt
#
# Or, if called from a main pipeline script:
#   sbatch 5_psiclass_assembly.sbatch "$GENOME_FASTA" "$REFERENCE_GTF" "$PSICLASS_OUTPUT_DIR" "$BAM_LIST_FILE"

# --- SLURM Directives ---
#SBATCH --job-name=psiclass_assembly     # Job name
#SBATCH --nodes=1                        # Request one node
#SBATCH --ntasks-per-node=1              # Request one task (process) per node
#SBATCH --cpus-per-task=16               # Request 16 CPU cores for this task (adjust as needed for PsiCLASS)
#SBATCH --mem=64G                        # Request 64 GB of memory (adjust based on genome size & sample count)
#SBATCH --time=12:00:00                  # Max job time (e.g., 12 hours, adjust as needed)
#SBATCH --output=slurm_psiclass_%j.out   # Standard output file
#SBATCH --error=slurm_psiclass_%j.err    # Standard error file
#SBATCH --partition=long                 # Specify a partition (e.g., 'long', 'medium', 'batch')
#SBATCH --mail-type=END,FAIL             # Email notifications for job end or failure
#SBATCH --mail-user=your.email@example.com # Your email address

# --- Load necessary modules (adjust based on your cluster's setup) ---
module purge
# module load gcc/9.3.0 # Example compiler module if PsiCLASS was compiled with it
# module load Python/3.11.5-GCCcore-13.2.0
# module load PsiCLASS/version_number # If PsiCLASS is available as a module
# If PsiCLASS is not a module, ensure it's in your PATH or provide full path
# export PATH="/path/to/psiclass/bin:$PATH" # Example if manually managing PATH
# module load SAMtools/1.19.2-GCC-13.2.0


WORK_DIR="/home/data/human_genetics/USERS/jhealy/rna_pipeline_v0.1"
REFERENCE_DIR="${WORK_DIR}/reference_data/human/GRCh38" 
GTF_FILE="${REFERENCE_DIR}/gencode.v48.chr_patch_hapl_scaff.annotation.gtf" 
# REFERENCE_FASTA="$1"  
# GTF_FILE="$2"    
# OUTPUT_DIR="$3"       
# BAM_LIST_FILE="$4"    

# --- Validate Inputs ---
if [ -z "$REFERENCE_FASTA" ] || [ -z "$GTF_FILE" ] || [ -z "$OUTPUT_DIR" ] || [ -z "$BAM_LIST_FILE" ]; then
    echo "Error: Missing required arguments." >&2
    echo "Usage: $0 <reference_fasta> <reference_gtf> <output_dir> <bam_list_file>" >&2
    exit 1
fi

if [ ! -f "$REFERENCE_FASTA" ]; then
    echo "Error: Reference FASTA file not found: $REFERENCE_FASTA" >&2
    exit 1
fi
if [ ! -f "$GTF_FILE" ]; then
    echo "Error: Reference GTF/GFF3 file not found: $GTF_FILE" >&2
    exit 1
fi
if [ ! -f "$BAM_LIST_FILE" ]; then
    echo "Error: BAM list file not found: $BAM_LIST_FILE" >&2
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR" || { echo "Error: Could not create output directory: $OUTPUT_DIR" >&2; exit 1; }

# --- Set PsiCLASS specific parameters ---
# Adjust these based on PsiCLASS documentation and your data
NUM_THREADS="$SLURM_CPUS_PER_TASK" # Use SLURM allocated CPUs
# For PsiCLASS, the --num-cpu argument is important.
# Check PsiCLASS documentation for other parameters:
# e.g., --min-isoform-len, --min-read-count, --min-exon-len, --max-intron-len
# PSI_CLASS_ARGS="--num-cpu ${NUM_THREADS} --min-isoform-len 100 --min-read-count 5" # Example additional arguments

# --- Execute PsiCLASS ---
echo "--- Starting PsiCLASS Assembly ---"
echo "Reference FASTA: $REFERENCE_FASTA"
echo "Reference GTF: $GTF_FILE"
echo "Output Directory: $OUTPUT_DIR"
echo "BAM List File: $BAM_LIST_FILE"
echo "Threads: $NUM_THREADS"
echo "Running command:"
echo "PsiCLASS -c $NUM_THREADS -g $GTF_FILE -f $REFERENCE_FASTA -o $OUTPUT_DIR -b $BAM_LIST_FILE"

# The core PsiCLASS command
PsiCLASS -c "${NUM_THREADS}" \
         -g "${GTF_FILE}" \
         -f "${REFERENCE_FASTA}" \
         -o "${OUTPUT_DIR}" \
         -b "${BAM_LIST_FILE}"

# Check exit status of PsiCLASS
if [ $? -eq 0 ]; then
    echo "--- PsiCLASS Assembly Completed Successfully! ---"
    echo "Output files are in: $OUTPUT_DIR"
else
    echo "--- Error: PsiCLASS Assembly Failed! ---" >&2
    exit 1
fi

echo "Job finished at: $(date)"

sbatch 5_psiclass_assembly.sbatch \
    <REFERENCE_FASTA_PATH> \
    <REFERENCE_GTF_PATH> \
    <PSICLASS_OUTPUT_DIR> \
    <BAM_LIST_FILE_PATH>