#!/bin/bash

#SBATCH --time=5:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=2G
#SBATCH --job-name="ExtractionPipeline"

# Here load needed system tools (Java 1.8 is strictly required, one of two: singularity or conda for python 2.7 are needed,
# depending on the method which is used for dependency management)
module load jdk/16.0.1
module load openjdk/11.0.2
module load squashfs
module load singularity

set -f

# NB! Following will extract all the genes and SNPs! In case of full mapping,
# it is advisable to use gene filter or P-value filter to get manageable subset!

nextflow_path=/gpfs/space/GI/eQTLGen/tools/

mastertable="/gpfs/space/GI/eQTLGen/freeze3/eqtl_mapping/input/mastertable_empirical_2024-08-29_extended.txt"

#empirical=/gpfs/space/GI/eQTLGen/freeze2/eqtl_mapping/output/empirical_4GenPCNoExpPC_2023-07-13_PerGene
empirical="/gpfs/space/GI/eQTLGen/freeze3/eqtl_mapping/output/all_empirical_4GenPCNoExpPC_2024-09-05/eqtls/"
genes="${empirical}/available_genes.txt"

echo "ID" > ${genes}
ls ${empirical}/meta | awk -F'=' '{print $2}' >> ${genes}

genome_reference="hg38.genome"
variant_reference="/gpfs/space/GI/eQTLGen/freeze1/eqtl_mapping/MetaAnalysis/bin/hase/data/1000G-30x.ref.gz"
gene_reference="/gpfs/space/GI/eQTLGen/freeze1/InputFilesForPaper/2023-01-28_MetaAnalysis/data/Homo_sapiens.GRCh38.106.gtf.gz"

inclusion_step_output="/gpfs/space/GI/eQTLGen/freeze3/InclusionLists/output_2024-08-29"

output_folder="/gpfs/space/GI/eQTLGen/freeze3/Interpretation/extractions/all_empirical_4GenPC20ExpPC_2024-09-05"

p_threshold="5e-8"

NXF_VER=21.10.6 ${nextflow_path}/nextflow run ../ExtractResults.nf \
--input ${empirical}/meta \
--genes ${genes} \
--variant_reference ${variant_reference} \
--gene_reference ${gene_reference} \
--mastertable ${mastertable} \
--p_threshold ${p_threshold} \
--output ${output_folder} \
--inclusion_step_output ${inclusion_step_output} \
-resume \
-profile slurm,singularity
