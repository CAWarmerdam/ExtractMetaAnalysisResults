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
Pipeline for parallelized extraction and filtering of the raw HASE results.

This pipeline is used to extract subsets of results from the HASE results (numerous large .parquet files).

Usage:

nextflow run ExtractHaseResults.nf \
--empirical '/inputfolder/' \
--permuted '/outputfile/' \
--genes '/phenotypes.txt' \
--ld-dataset '/dataset/' \
--output '/output/'


Mandatory arguments:
--empirical           Path to the folder with HASE result .parquet files.
--permuted            Path to where the database should be written
--genes               Path to a file with all unique genes
--maf-table           Path to table with maf per variant
--ld-dataset          Path to LD dataset
--output              Path to outputfolder

""".stripIndent()

}

if (params.help){
    helpmessage()
    exit 0
}


// Define list of chromosomes to analyse
chromosomes = params.chromosome ? [params.chromosome.toString()] : (1..22).collect { it.toString() }
//chromosomes = [1,22]

log.info(chromosomes)

//Default parameters
Channel.fromPath(params.dataset).collect().set { dataset_parquet_ch }
Channel.fromPath(params.pca_folder).collect().set { pca_folder_ch }
Channel.fromPath(params.genes).collect().set { ld_genes_ch }
Channel.fromPath(params.variant_reference).collect().set { variant_reference_ch }
Channel.fromList(chromosomes).map( chr -> tuple(chr)).view().set { chromosome_ch }

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
summary['Empirical eQTLs']                          = params.empirical
summary['LD data']                                  = params.ld
summary['Significant eQTLs']                        = params.significant
summary['Genome reference']                         = params.genome_reference
summary['Variant reference']                        = params.variant_reference
summary['Gene reference']                           = params.gene_reference
summary['Gene list']                                = params.genes
summary['Chromosomes']				    = chromosomes

log.info summary.collect { k,v -> "${k.padRight(21)}: $v" }.join("\n")
log.info "======================================================="

workflow {
    // Buffer genes
    variant_chunk_ch = SplitVariantSet(variant_reference_ch, number_of_chunks).splitCsv( header: false, skip: 1 ).combine(chromosome_ch, by: 0).view()
    GenerateLdPanel(dataset_parquet_ch, pca_folder_ch, variant_reference_ch, ld_genes_ch, variant_chunk_ch)
}

workflow.onComplete {
    println ( workflow.success ? "Pipeline finished!" : "Something crashed...debug!" )
}
