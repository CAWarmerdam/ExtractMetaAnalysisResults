#!/bin/bash

#SBATCH --time=1-00:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=2G
#SBATCH --job-name="LdPipeline"

# Here load needed system tools (Java 1.8 is strictly required, one of two: singularity or conda for python 2.7 are needed,
# depending on the method which is used for dependency management)
module load jdk/16.0.1
module load openjdk/11.0.2
module load squashfs
module load singularity

set -f

# NB! Following will extract all the genes and SNPs! In case of full mapping,
# it is advisable to use gene filter or P-value filter to get manageable subset!

nextflow_path=/gpfs/space/GI/eQTLGen/freeze2/tools

empirical=/gpfs/space/GI/eQTLGen/freeze2/eqtl_mapping/output/empirical_4GenPC20ExpPC_2023-05-10_perGene/
permuted=/gpfs/space/GI/eQTLGen/freeze1/eqtl_mapping/output/permuted_4GenPC20ExpPC_2023-04-18_PerGene/

reference_data='/gpfs/space/GI/eQTLGen/make_reference_files/30x_reference/hg38/phasing_reference/phasing/chr*.bcf'

genome_reference="/gpfs/space/GI/eQTLGen/freeze1/Interpretation/ld/ExtractMetaAnalysisResults/hg38.genome"
variant_reference="/gpfs/space/GI/eQTLGen/make_reference_files/30x_reference/hg38/hase_reference/1000G-30x.ref.gz"
gene_reference="/gpfs/space/GI/eQTLGen/freeze1/Interpretation/ld/ExtractMetaAnalysisResults/gencode.v43.basic.annotation.gff3.gz"

bed="/gpfs/space/GI/eQTLGen/freeze1/Interpretation/ld/data/vuckovic_flanked1mb_na_significant_loci_5mlog9_hg38.bed"

inclusion_step_output="/gpfs/space/GI/eQTLGen/freeze2/InclusionLists/output"
maf_table="/gpfs/space/GI/eQTLGen/freeze2/PreMetaQc/output/MafInformation.txt.gz"

output_folder="/gpfs/space/GI/eQTLGen/freeze2/Interpretation/ld/output2"

NXF_VER=23.04.1 ${nextflow_path}/nextflow run ExtractMainResults.nf \
--empirical ${empirical}/eqtls \
--permuted ${permuted}/eqtls \
--reference_data ${reference_data} \
--genes ${empirical}/phenotypes_unique.txt \
--genome_reference ${genome_reference} \
--variant_reference ${variant_reference} \
--gene_reference ${gene_reference} \
--background_bed ${bed} \
--maf_table ${maf_table} \
--output ${output_folder} \
--inclusion_step_output ${inclusion_step_output} \
-resume \
-profile slurm,singularity
