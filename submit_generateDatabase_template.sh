#!/bin/bash

#SBATCH --time=1:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=2G
#SBATCH --job-name="GenerateDatabase"

# Here load needed system tools (Java 1.8 is strictly required, one of two: singularity or conda for python 2.7 are needed,
# depending on the method which is used for dependency management)
module load java/11.0.8_10
module load any/singularity/3.5.3
module load squashfs/4.4

# NB! Following will extract all the genes and SNPs! In case of full mapping,
# it is advisable to use gene filter or P-value filter to get manageable subset!

nextflow_path=/gpfs/space/GI/eQTLGen/EstBB_testing/MetaAnalysis/sqliteTest/2022-11-17/tools/

input_folder=/gpfs/space/GI/eQTLGen/freeze1/eqtl_mapping/output/empirical_4GenPC20ExpPC_2022-11-14/MetaAnalysisResultsEncoded/
output_file=eqtlgen_metaanalysis_2022_11_17.db
db_folder=/gpfs/space/GI/eQTLGen/EstBB_testing/MetaAnalysis/sqliteTest/2022-11-17/db
output_folder=/gpfs/space/GI/eQTLGen/EstBB_testing/MetaAnalysis/sqliteTest/2022-11-17/output

NXF_VER=20.10.0 ${nextflow_path}/nextflow run generateDatabase.nf \
--inputfolder ${input_folder} \
--parquet 'node_27_*_??_result.parquet' \
--dbfile ${output_file} \
--dbfolder ${db_folder} \
--outputfolder ${output_folder} \
-resume \
-profile slurm,singularity
