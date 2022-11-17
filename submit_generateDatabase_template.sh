#!/bin/bash

#SBATCH --time=24:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=2G
#SBATCH --job-name="ExtractHaseResults"

# Here load needed system tools (Java 1.8 is strictly required, one of two: singularity or conda for python 2.7 are needed,
# depending on the method which is used for dependency management)
module load any/jdk/1.8.0_265
module load any/singularity/3.5.3
module load squashfs/4.4

# NB! Following will extract all the genes and SNPs! In case of full mapping, 
# it is advisable to use gene filter or P-value filter to get manageable subset!

nextflow_path=[path to folder where Nextflow executable is]

input_folder=[folder with .parquet files]
output_file=[output database file]

NXF_VER=20.10.0 ${nextflow_path}/nextflow run generateDatabase.nf \
--inputfolder ${input_folder} \
--outputfile ${output_file} \
-resume \
-profile slurm,singularity