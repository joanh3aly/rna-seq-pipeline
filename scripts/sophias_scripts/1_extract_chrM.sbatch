#!/bin/bash
####
#### Sophia Heneghan <sophiaheneghan22@rcsi.ie>
#### hg38_Mitochondrial_Variant_calling
####
#SBATCH --job-name=chrM_extract
#SBATCH --output=/home/data/human_genetics/SHeneghan/Mitochondrial_Variants/slurmfiles/%j.chrM_extract.out
#SBATCH --error=/home/data/human_genetics/SHeneghan/Mitochondrial_Variants/slurmfiles/%j.chrM_extract.err
#SBATCH --partition=ProdQ
#SBATCH --time=2-00:00:00
#SBATCH --ntasks=1
#SBATCH --mem=4000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=sophiaheneghan22@rcsi.com
#SBATCH --array=0-277%10

WORK_DIR=/home/data/human_genetics/SHeneghan/Mitochondrial_Variants
REF_GENOME=/home/data/lighthouse/suppfiles/hg38_v0_Homo_sapiens_assembly38.fasta
EXOME_BAM=/home/data/human_genetics/Renal/Exomes/CeGaT_PKD/BAM
TGP_BAM=/home/data/human_genetics/Renal/Joint_Calling/BAM
OUT_DIR=$WORK_DIR/1_chrM_extracted
EXOME_REGIONS=/home/data/human_genetics/reference/hg38/S000021_hg38_targets_extended.bed
TGP_REGIONS=/home/data/human_genetics/reference/hg38/171113_HG38_Renal227_EZ_HX3_capture_targets.bed

###go through directory and extract filepaths for each bam file

module load GATK/4.2.0.0-GCCcore-10.2.0-Java-1.8.0_202

for f in "$EXOME_BAM"/*; do
    if [[ $f =~ ^${EXOME_BAM}/.*.markdup.recal.bam$ ]]; then
      arrVar+=("$f")
    fi
done


FILE="${arrVar[$SLURM_ARRAY_TASK_ID]}"
###FILE is the one that prints each of the BAM files

echo $FILE

tmp=${FILE#*BAM/}
SAMPLE=${tmp%.markdup.recal.bam*}

echo $SAMPLE

gatk PrintReads \
	-I $FILE \
	-L chrM \
	-O ${WORK_DIR}/1_chrM_extracted/${SAMPLE}_chrM.bam

