#!/bin/bash

#SBATCH --time=1-12:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=4G
#SBATCH --job-name="LdPipeline"

# Here load needed system tools (Java 1.8 is strictly required, one of two: singularity or conda for python 2.7 are needed,
# depending on the method which is used for dependency management)
#module load jdk/16.0.1
#module load openjdk/11.0.2
#module load squashfs
#module load singularity
module load Java/11.0.16 # update java version

set -f

# NB! Following will extract all the genes and SNPs! In case of full mapping,
# it is advisable to use gene filter or P-value filter to get manageable subset!

nextflow_path=/scratch/hb-functionalgenomics/projects/eqtlgen-phase2/nextflow

empirical=/scratch/hb-functionalgenomics/projects/eqtlgen-phase2/freeze3/eqtl_mapping/output/2024-09-30_meta_analysis/
#permuted=/scratch/hb-functionalgenomics/projects/eqtlgen-phase2/freeze3/eqtl_mapping/output/2024-10-09_permuted/
permuted="/scratch/hb-functionalgenomics/projects/eqtlgen-phase2/fine_mapping/jobs/process_ld_panel/ld_panel"

#significant="/scratch/hb-functionalgenomics/projects/eqtlgen-phase2/freeze3/Interpretation/cis-trans-coloc/input/subset_p5e8_hyprColocFormat_2024-09-05.csv.gz"
significant="/scratch/hb-functionalgenomics/projects/eqtlgen-phase2/freeze3/Interpretation/extractions/all_empirical_4GenPC20ExpPC_2024-09-05/loci_empirical_annotated/2024-09-05_all_genes_sign_5em8_cis_window_1.5e6_annotated.hyprcoloc_format.txt.gz"
#significant="/scratch/hb-functionalgenomics/projects/eqtlgen-phase2/freeze3/Interpretation/extractions/all_empirical_4GenPC20ExpPC_2024-09-05/loci_empirical_annotated/2024-09-05_all_genes_sign_5em8_cis_window_1.5e6_annotated_chr1.hyprcoloc_format.txt.gz"

uncorrelated_genes="/scratch/hb-functionalgenomics/projects/eqtlgen-phase2/fine_mapping/input/uncorrelated_gene_list_variable_availability_n2051_rsq0.02_genes_available_20241127.txt"

genome_reference="/scratch/hb-functionalgenomics/projects/eqtlgen-phase2/public_data/hg38.genome"
variant_reference="/scratch/hb-functionalgenomics/projects/eqtlgen-phase2/processed_data/variants/1000G-30x_index.parquet"
gene_reference="/scratch/hb-functionalgenomics/projects/eqtlgen-phase2/public_data/Homo_sapiens.GRCh38.106.gtf.gz"

output_folder="/scratch/hb-functionalgenomics/projects/eqtlgen-phase2/fine_mapping/output_alt"

# Use code below to generate a list of all genes:
# Make sure this is a file with a header, and genes (ENSG...). Only available genes, only those will be tested. If more genes, an error will be thrown
#echo "ID" > ${empirical}/available_genes.txt
ls ${empirical}/meta | awk -F'=' '{print $2}' > ${empirical}/available_genes.txt

# Use code below to generate a list of all uncorrelated genes:
# Make sure this is a file with a header, and genes (ENSG...). Only available genes, only those will be tested. If more genes, an error will be thrown
#echo "ID" > ${permuted}/uncorrelated_genes.txt
#ls ${permuted}/eqtls/meta | awk -F'=' '{print $2}' > ${permuted}/uncorrelated_genes.txt

NXF_VER=23.04.1 ${nextflow_path}/nextflow run FineMapping.nf \
 --empirical ${empirical}/meta \
 --permuted ${permuted} \
 --significant ${significant} \
 --genes ${empirical}/available_genes.txt \
 --uncorrelated_genes ${uncorrelated_genes} \
 --genome_reference ${genome_reference} \
 --variant_reference ${variant_reference} \
 --gene_reference ${gene_reference} \
 --output ${output_folder} \
 --chromosome 21 \
 -resume \
 -work-dir work_test \
 -profile slurm,singularity
