#!/bin/bash
####
#### Joan Healy <joanhealy24@rcsi.ie>
#### Index reference genome
#### Inputs: ref genome and gtf 
#### Outputs: ref genome index files
####
#SBATCH --job-name=hisat2_index_GRCh38
#SBATCH --output=/home/data/human_genetics/USERS/jhealy/rna_pipeline_v0.1/reference_data/human/GRCh38/HISAT2_index/slurmfiles/%j_pipeline_v0.1.out
#SBATCH --error=/home/data/human_genetics/USERS/jhealy/rna_pipeline_v0.1/reference_data/human/GRCh38/HISAT2_index/slurmfiles/%j_pipeline_v0.1.err
#SBATCH --partition=ProdQ 
#SBATCH --time=7-12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16 
#SBATCH --mem=350G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joanhealy24@rcsi.com
# #SBATCH --nodelist=prdubrhpc1-05
# #SBATCH --reservation=TAPM


set -x
# Load modules
module load HISAT2/2.2.1-foss-2019b
module load Python/3.11.5-GCCcore-13.2.0

# Define paths
REFERENCE_DIR="/home/data/human_genetics/USERS/jhealy/rna_pipeline_v0.1/reference_data/human/GRCh38" 

# HISAT2_INDEX_DIR="${REFERENCE_DIR}/HISAT2_index/v_2.2.1-foss-2019b_basic"
# HISAT2_INDEX_DIR="${REFERENCE_DIR}/HISAT2_index/p14_masked"
HISAT2_INDEX_DIR="${REFERENCE_DIR}/HISAT2_index/p14_masked_gtf_masked"

# REFERENCE_FASTA="${REFERENCE_DIR}/GRCh38.p14.genome.fa"
REFERENCE_FASTA="${REFERENCE_DIR}/masked_basic_genome.fasta"
# REFERENCE_FASTA="${REFERENCE_DIR}/masked_genome.fasta"
# REFERENCE_FASTA="${REFERENCE_DIR}/masked_pc_trans_genome_2.fasta"

# GTF_FILE="${REFERENCE_DIR}/gencode.v48.chr_patch_hapl_scaff.annotation.gtf" 
# GTF_FILE="${REFERENCE_DIR}/gencode.v48.basic.annotation_no_pseudogenes.gtf" 
GTF_FILE="${REFERENCE_DIR}/gencode.v48.basic.annotation_no_PKDPs_2.gtf"

SNP_VCF_GZ_ORIGINAL="${REFERENCE_DIR}/00-common_all.vcf.gz"
SNP_VCF_FILE="${REFERENCE_DIR}/00-common_all.vcf"
SNP_VCF_CHRFX_FILE="${REFERENCE_DIR}/00-common_all_chr_fixed.vcf" 


# Directory where the index files will be created
mkdir "${HISAT2_INDEX_DIR}"
cd "${HISAT2_INDEX_DIR}"

printf "Starting HISAT2 index building for GRCh38...\n"
printf "Reference DIR: %s\n" "${REFERENCE_DIR}"
printf "Reference FASTA: %s\n" "${REFERENCE_FASTA}"
printf "GTF File: %s\n" "${GTF_FILE}"
printf "SNP zip File: %s\n" "${SNP_VCF_GZ_ORIGINAL}"
printf "SNP File: %s\n" "${SNP_VCF_FILE}"
printf "Output Index Directory: %s\n" "${HISAT2_INDEX_DIR}"
printf "Chr-fixed SNP File: %s\n" "${SNP_VCF_CHRFX_FILE}"



# Extract splice sites 
printf "Extracting splice sites...\n"

HISAT2_SCRIPTS_DIR="${EBROOTHISAT2}/bin"
if [ ! -d "$HISAT2_SCRIPTS_DIR" ]; then
    printf "WARNING: Could not find HISAT2 scripts directory. Please check if HISAT2 was installed in a standard way or manually specify HISAT2_SCRIPTS_DIR.\n"
fi

if [ -f "${HISAT2_SCRIPTS_DIR}/extract_splice_sites.py" ]; then
    python3 "${HISAT2_SCRIPTS_DIR}/extract_splice_sites.py" "${GTF_FILE}" > genome.ss
    if [ $? -ne 0 ]; then echo "Error extracting splice sites. Exiting." >&2; exit 1; fi
else
    echo "ERROR: extract_splice_sites.py not found at ${HISAT2_SCRIPTS_DIR}. Cannot build splice-aware index. Exiting." >&2
    exit 1
fi

printf "Extracting exons...\n"
if [ -f "${HISAT2_SCRIPTS_DIR}/extract_exons.py" ]; then
    python3 "${HISAT2_SCRIPTS_DIR}/extract_exons.py" "${GTF_FILE}" > genome.exon
    if [ $? -ne 0 ]; then echo "Error extracting exons. Exiting." >&2; exit 1; fi
else
    echo "ERROR: extract_exons.py not found at ${HISAT2_SCRIPTS_DIR}. Cannot build splice-aware index. Exiting." >&2
    exit 1
fi

# printf "Extracting SNPs...\n"
# if [ -f "${HISAT2_SCRIPTS_DIR}/hisat2_extract_snps_haplotypes_VCF.py" ]; then
#     # zcat to decompress and pipe to awk for chr prefixing
#     zcat "${SNP_VCF_GZ_ORIGINAL}" | awk 'BEGIN {OFS="\t"} {if ($0 !~ /^#/) {$1="chr"$1} print $0}' > "${SNP_VCF_CHRFX_FILE}"

#     if [ $? -ne 0 ]; then
#         echo "ERROR: Failed to fix VCF chromosome names. Exiting." >&2
#         exit 1
#     fi
#     # Ensure the newly created file exists
#     if [ ! -f "${SNP_VCF_CHRFX_FILE}" ]; then
#         echo "ERROR: Chr-fixed VCF file ${SNP_VCF_CHRFX_FILE} was not created. Exiting." >&2
#         exit 1
#     fi
#     printf "VCF chromosome names fixed: ${SNP_VCF_CHRFX_FILE}\n"

#     # Use the chr-fixed VCF file directly
#     # python3 "${HISAT2_SCRIPTS_DIR}/hisat2_extract_snps_haplotypes_VCF.py" "${SNP_VCF_CHRFX_FILE}" > genome.snp
#     python3 "${HISAT2_SCRIPTS_DIR}/hisat2_extract_snps_haplotypes_VCF.py" \
#         "${REFERENCE_FASTA}" \
#         "${SNP_VCF_CHRFX_FILE}" \
#         genome
#     if [ $? -ne 0 ]; then
#         echo "ERROR: Failed to extract SNPs from chr-fixed VCF ${SNP_VCF_CHRFX_FILE}. Exiting." >&2
#         exit 1
#     fi
# else
#     echo "ERROR: hisat2_extract_snps_haplotypes_VCF.py not found at ${HISAT2_SCRIPTS_DIR}. Cannot build SNP-aware index. Exiting." >&2
#     exit 1
# fi
# printf "SNP extraction process finished (output to genome.snp).\n"

printf "Building HISAT2 index with splice sites and exons...\n"
printf "Value of SLURM_CPUS_PER_TASK: '%s'\n" "${SLURM_CPUS_PER_TASK}"
if [ -z "${SLURM_CPUS_PER_TASK}" ]; then
    NUM_THREADS=16
    printf "WARNING: SLURM_CPUS_PER_TASK is empty. Using default %d threads for hisat2-build.\n" "${NUM_THREADS}"
else
    NUM_THREADS=${SLURM_CPUS_PER_TASK}
    printf "Using %d threads for hisat2-build as specified by SLURM_CPUS_PER_TASK.\n" "${NUM_THREADS}"
fi

# hisat2-build -p ${NUM_THREADS} --ss ss.txt --exon exon.txt "${REFERENCE_FASTA}" GRCh38_splices_index
# hisat2-build -p ${NUM_THREADS}  \
#             --snp genome.snp \ 
#             --haplotype genome.haplotype \
#             --exon genome.exon \
#             --ss genome.ss \
#             "${REFERENCE_FASTA}" \
#             GRCh38_splices_index_snp_tran

# hisat2-build -p ${NUM_THREADS}  \
#             --exon genome.exon \
#             --ss genome.ss \
#             "${REFERENCE_FASTA}" \
#             GRCh38_splices_index_snp_tran

# hisat2-build -p 8  \
#             --snp genome.snp \ 
#             --exon genome.exon \
#             --ss genome.ss \
#             "${REFERENCE_FASTA}" \
#             GRCh38_splices_index_snp_tran

module load HISAT2/2.2.1-foss-2019b
hisat2-build -p 16 \
            --exon genome.exon \
            --ss genome.ss \
            "${REFERENCE_FASTA}" \
            GRCh38_basic_ss_exon_index
            # --snp genome.snp \
            # --haplotype genome.haplotype \
            #  --large-index 


if [ $? -eq 0 ]; then
    printf "\nHISAT2 index building completed successfully for GRCh38_splices_index.\n"
else
    printf "\nERROR: HISAT2 index building failed.\n" >&2
    exit 1
fi

# Clean up temporary files
# rm ss.txt exon.txt

printf "Finished HISAT2 index building script.\n"









