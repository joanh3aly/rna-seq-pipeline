# config.sh
REFERENCE_DIR="/home/data/human_genetics/USERS/jhealy/rna_pipeline_v0.1/reference_data"
GENOME_FASTA="${REFERENCE_DIR}/human/GRCh38/GRCh38.p14.genome.fa"
REFERENCE_GTF="${REFERENCE_DIR}/annotations/gencode.v44.annotation.gtf"

# Base output directory for this pipeline run
PIPELINE_OUTPUT_ROOT="/scratch/jhealy/my_rna_seq_project_$(date +%Y%m%d_%H%M%S)" # Unique run dir

# Define stage-specific parent directories (will be created by master script or relevant stage script)
RAW_READS_DIR="${PIPELINE_OUTPUT_ROOT}/01_raw_reads"
TRIMMED_READS_DIR="${PIPELINE_OUTPUT_ROOT}/02_trimmed_reads"
ALIGNED_BAMS_DIR="${PIPELINE_OUTPUT_ROOT}/03_aligned_bams"
TRANSCRIPT_ASSEMBLY_DIR="${PIPELINE_OUTPUT_ROOT}/04_transcript_assembly"
# ... and so on