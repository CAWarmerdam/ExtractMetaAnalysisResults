#!/bin/bash

#SBATCH --time=1:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=2G
#SBATCH --job-name="GenerateDatabase"

# Here load needed system tools (Java 1.8 is strictly required, one of two: singularity or conda for python 2.7 are needed,
# depending on the method which is used for dependency management)
module load jdk/16.0.1
module load openjdk/11.0.2
module load any/singularity
module load squashfs/4.4

# NB! Following will extract all the genes and SNPs! In case of full mapping,
# it is advisable to use gene filter or P-value filter to get manageable subset!

nextflow_path=/gpfs/space/GI/eQTLGen/EstBB_testing/MetaAnalysis/tools

input_folder=/gpfs/space/GI/eQTLGen/hase_output_testing/output/freeze1/eqtls/
output_folder=/gpfs/space/GI/eQTLGen/hase_output_testing/output/t
variants=/gpfs/space/GI/eQTLGen/hase_output_testing/sampled_variants.txt

NXF_VER=21.10.6 ${nextflow_path}/nextflow run analysis.nf \
--inputfolder ${input_folder} \
--outputfolder ${output_folder} \
--variants ${variants} \
-resume \
-profile slurm,singularity
