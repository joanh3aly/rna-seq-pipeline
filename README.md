# Splice site Transcriptomics Pipeline
This repo is a transcriptomics pipeline for RNA seq data, created for the Futureneuro lab in the Royal College of Surgeons Ireland.
It was specifically designed to use an index customised for the discovery of splice site variants in paired end RNA-seq data. The pipeline was created using Slurm sbatch scripts.

## There are 6 steps in the pipeline:
1. Creation of run ID, run directory, all subdirectories, and the copying of fastq files into the correct subdirectory.
2. Merging of all sequence sample lanes into one FASTQ file
3. FastQC on raw FASTQ files
4. Trimming and filtering of low quality FASTQ files using Trimmomatic
5. FastQC of trimmed files to check for improved quality
6. Alignment with the HISAT2 index


## Scripts:
All sbatch scripts are located in the scripts directory.
- Entry script - pipeline_v0.1.sbatch
- Parallelisation script (for faster processing) - parallelize_fastq_module.sbatch
- Directory and run ID generation script - 
- Concatenate sequence lanes script - 
- QC raw fastq script - 1_fastqc_untrimmed.sbatch
- Trim script - 2_trim_fastq_trimmomatic.sbatch or 2_trim_fastq_bbmap.sbatch
- QC trimmed fastq - 3_fastqc_trimmed.sbatch
- Alignment script - 4.hisat2_align.sbatch
- To split fastq files that are interleaved - split_interleaved.sh


## How to run:
1. Add all fastq files into the "input_drop" directory.  
IMPORTANT: Fastq files must have the correct format. The format we used is  `<rnaseqID>_<lane>_<read>.fq.gz`. For example: SRR10488341_L1_1.fq.gz (lane 1, read 1) /SRR10488341_L1_2.fq.gz (lane 1, read 2).  
There must be a paired end file for the scripts in the pipeline to work.  
Both .fq and .fastq extensions are handled by the regex in the script.  
2. Add the correct working directory filepath to the top of the pipeline_v0.1 directory for the WORK_DIR constant.
3. Run the entry script in the bash console: sbatch pipeline_v0.1.sbatch  
This will run each of the 6 stages of the pipeline in the correct order.  
The run ID directory will be automatically created in the base directory from the current date and time of the run  `run_<yearmonthdate>_<hourminute>` - ie run_20250618_2125


## Indexing:
- The HISAT2 indexing script runs independently from the pipeline. 
- hisat2_index_ref.sbatch can be used to create a new index. It is located in the scripts directory.
- Run sbatch hisat2_index_ref.sbatch in the bash console to execute this script.
- All indices, reference genomes and supporting files (gtf, bed, vcf etc) are in the reference_data folder.
- All output files from the indexing script should be saved in the  `reference_data/human/GRCh38/HISAT2_index/v_2.2.1-foss-2019b_<name of custom version>` directory.
    * The indexing script should create 8 ht2 files (these are used by the alignment step) and 4 temporary files that are used to generate the ht2 files (genome.exon/genome.haplotype/genome.snp/genome.ss)
- Currently there are 2 human genome index options in the human/GRCH38 subdirectories:
    * reference_data/human/GRCh38/HISAT2_index/v_2.2.1-foss-2019b_pc_trans_masked - this index was generated from the exon and ss (splice site) temp files, with blacklisted PKD1 pseudogenes from the protein coding GRCH38 reference genome (masked_pc_trans_genome_2.fasta) and basic gtf file (gencode.v48.basic.annotation_no_pseudogenes.gtf).
        * Masking/blacklisting of pseudogenes was completed using BEDtools.
    * reference_data/human/GRCh38/HISAT2_index/v_2.2.1-foss-2019b contains an index created from ss and exons.
    * See the gencode website for more reference genome and gtf options: https://www.gencodegenes.org/human/
- A prebuilt index from the HISAT2 website is available in the grch38_snp directory. This contains SNPs, haplotypes, exons, and splice site data. https://daehwankimlab.github.io/hisat2/download/#h-sapiens


## Log files for debugging:
- All error and output files for each step in the pipeline run are put into the  `<runIDdirectory>/slurmfiles` directory in each unique run directory.
- For the pipeline entry script itself, the err/out files are saved in scripts/slurmfiles. This is because there can be pipeline script errors before the run's directories are created.
- Err and out files for the indexing script are saved into rna_pipeline_v0.1/reference_data/human/GRCh38/HISAT2_index/slurmfiles, as this is a separate process from the pipeline.
- Trimmomatic creates logs - these are also in the  `<runIDdirectory>/slurmfiles` directory.
- HISAT2 creates logs for the alignment step - these files are in the aligner_outputs directory.


## Improvements:
- Automated unit tests
- New options for choosing different aligners ie STAR
- Options to create a new custom index as part of the pipeline.
- Make indexing more user friendly (input options for genome version, gtf options etc)
- Assembly step - ie Psiclass script was in progress but never finished
- Visualisation step - ie: using Spladder for mathemathical splicing analysis. https://github.com/ratschlab/spladder
- Splice site protein prediction - Bisbee https://github.com/tgen/bisbee
- Making the processing of input files more flexible/robust for different FASTQ file naming conventions
- Better error handling
- IaC, CI/CD pipelines for automated environment deployments.


## Versions of all modules used:
FastQC/0.11.9-Java-11  
BBMap/38.90-GCC-10.2.0  
R/4.4.1-gfbf-2023b  
Java/11.0.16  
Trimmomatic/0.39-Java-11  
Python/3.11.5-GCCcore-13.2.0  
HISAT2/2.2.1-foss-2019b  
SAMtools/1.19.2-GCC-13.2.0  
