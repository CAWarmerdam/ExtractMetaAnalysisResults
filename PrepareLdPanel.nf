#!/usr/bin/env nextflow

/*
 * enables modules
 */
nextflow.enable.dsl = 2

// import modules
include { SplitVariantSet; GenerateLdPanel } from './modules/LdPanel'

def helpmessage() {

log.info"""

HASE output analyzer v${workflow.manifest.version}
==============================================
Pipeline for making an LD panel from permuted summary statistics and principal component analysis.
This used the PCA analyis performed in 'scripts/gene_coexpression.R'. That in turn uses the set of uncorrelated variants
from ExtractUncorrelatedVariants.nf, as well as the gene-gene correlations from ExtractGeneCorrelations.nf

Usage:

nextflow run PrepareLdPanel.nf \
--dataset '/inputfolder/' \
--pca_folder_prefix '/outputfile/' \
--genes '/phenotypes.txt' \
--variant_reference '/dataset/' \
--output 'output_path'

Mandatory arguments:
--dataset             Path to the folder with permuted HASE result .parquet files.
--pca_folder_prefix   Path to where the PCA results are stored
--genes               Path to a file with all unique genes
--variant_reference   Path to the variant reference parquet file

""".stripIndent()

}

if (params.help){
    helpmessage()
    exit 0
}


// Define list of chromosomes to analyse
chromosomes = params.chromosome ? params.chromosome.toString().tokenize(',').collect { it.toString() } : (1..22).collect { it.toString() }
//chromosomes = [1,22]

log.info(chromosomes.join(" "))

def pca_folder_prefix_file_object = file(params.pca_folder_prefix)
def pca_folder = pca_folder_prefix_file_object.parent
def pca_prefix = pca_folder_prefix_file_object.name

//Default parameters
Channel.fromPath(params.dataset).collect().set { dataset_parquet_ch }
Channel.fromPath(pca_folder).collect().set { pca_folder_ch }
Channel.fromPath(params.genes).collect().set { ld_genes_ch }
Channel.fromPath(params.variant_reference).collect().set { variant_reference_ch }
Channel.fromList(chromosomes).map( chr -> tuple(chr)).set { chromosome_ch }

number_of_chunks=200

log.info """=======================================================
HASE output analyzer v${workflow.manifest.version}"
======================================================="""
def summary = [:]
summary['Pipeline Version']                         = workflow.manifest.version
summary['Current user']                             = "$USER"
summary['Current home']                             = "$HOME"
summary['Current path']                             = "$PWD"
summary['Working dir']                              = workflow.workDir
summary['Script dir']                               = workflow.projectDir
summary['Config Profile']                           = workflow.profile
summary['Container Engine']                         = workflow.containerEngine
if(workflow.containerEngine) summary['Container']   = workflow.container
summary['Parquet dataset']                          = params.dataset
summary['PCA folder']                               = params.pca_folder_prefix
summary['Variant reference']                        = params.variant_reference
summary['Gene list']                                = params.genes
summary['Chromosomes']				                = chromosomes
summary['Output']				                    = params.output

log.info summary.collect { k,v -> "${k.padRight(21)}: $v" }.join("\n")
log.info "======================================================="

workflow {
    // Buffer genes
    variant_chunk_ch = SplitVariantSet(variant_reference_ch, number_of_chunks).splitCsv( header: false, skip: 1 ).combine(chromosome_ch, by: 0)
    GenerateLdPanel(dataset_parquet_ch, pca_folder_ch, variant_reference_ch, ld_genes_ch, variant_chunk_ch, pca_prefix)
}

workflow.onComplete {
    println ( workflow.success ? "Pipeline finished!" : "Something crashed...debug!" )
}
