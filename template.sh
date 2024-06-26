#!/bin/bash

#SBATCH --time=1-00:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=2G
#SBATCH --job-name="LdPipeline"

# Here load needed system tools (Java 1.8 is strictly required, one of two: singularity or conda for python 2.7 are needed,
# depending on the method which is used for dependency management)
#module load jdk/16.0.1
#module load openjdk/11.0.2
#module load squashfs
#module load singularity
#module load AdoptOpenJDK/11.0.4_11-hotspot
#module load nextflow/21.10.6-Java-11-LTS
module load Java/17.0.6

set -f

# NB! Following will extract all the genes and SNPs! In case of full mapping,
# it is advisable to use gene filter or P-value filter to get manageable subset!

nextflow_path=/scratch/hb-functionalgenomics/projects/eqtlgen-phase2/nextflow

empirical=/scratch/hb-functionalgenomics/projects/eqtlgen-phase2/output/2024-01-11-lld-GTEx-tmp-freeze/empirical_4GenPC20ExpPC_2024-02-25/
permuted=/scratch/hb-functionalgenomics/projects/eqtlgen-phase2/output/2024-01-11-lld-GTEx-tmp-freeze/permuted_4GenPC20ExpPC_2024-02-23/

reference_data="/scratch/hb-functionalgenomics/projects/eqtlgen-phase2/fine_mapping/ExtractMetaAnalysisResults/hg38/phasing_reference/phasing/chr*.bcf"

genome_reference="/scratch/hb-functionalgenomics/projects/eqtlgen-phase2/output/2024-01-11-lld-GTEx-tmp-freeze/ld_analysis_input/hg38.genome"
variant_reference="/scratch/hb-functionalgenomics/projects/eqtlgen-phase2/fine_mapping/ExtractMetaAnalysisResults/hg38/hase_reference/1000G-30x.ref.gz"
gene_reference="/scratch/hb-functionalgenomics/projects/eqtlgen-phase2/fine_mapping/ExtractMetaAnalysisResults/gencode.v43.basic.annotation.gff3.gz"

bed="/scratch/hb-functionalgenomics/projects/eqtlgen-phase2/output/2024-01-11-lld-GTEx-tmp-freeze/ld_analysis_input/vuckovic_flanked1mb_na_significant_loci_5mlog9_hg38.bed"

inclusion_step_output="/scratch/hb-functionalgenomics/projects/eqtlgen-phase2/output/2024-01-11-lld-GTEx-tmp-freeze/ld_analysis_input/freeze2/InclusionLists/output"
maf_table="/scratch/hb-functionalgenomics/projects/eqtlgen-phase2/output/2024-01-11-lld-GTEx-tmp-freeze/ld_analysis_input/freeze2/PreMetaQc/output/MafInformation_GTEx_LL.txt.gz"

#mastertable="/scratch/hb-functionalgenomics/projects/eqtlgen-phase2/output/2024-01-11-lld-GTEx-tmp-freeze/empirical_4GenPC20ExpPC_2023-10-30/MasterTable_2023-12-05.txt"
mastertable="/scratch/hb-functionalgenomics/projects/eqtlgen-phase2/output/2024-01-11-lld-GTEx-tmp-freeze/permuted_4GenPC20ExpPC_2023-10-30/MasterTable_permuted_2024-01-09.txt"

output_folder="/scratch/hb-functionalgenomics/projects/eqtlgen-phase2/fine_mapping/output"

#NXF_VER=23.04.1 ${nextflow_path}/nextflow run ExtractMainResults.nf \
${nextflow_path}/nextflow run FineMapping.nf \
--empirical ${empirical}/eqtls/meta \
--permuted ${permuted}/eqtls/meta \
--reference_data ${reference_data} \
--genes ${empirical}/phenotypes_unique.txt \
--genome_reference ${genome_reference} \
--variant_reference ${variant_reference} \
--gene_reference ${gene_reference} \
--background_bed ${bed} \
--maf_table ${maf_table} \
--output ${output_folder} \
--inclusion_step_output ${inclusion_step_output} \
--mastertable ${mastertable} \
-resume \
-profile slurm,singularity
#-profile singularity
