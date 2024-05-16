#!/usr/bin/env nextflow

/*
 * enables modules
 */
nextflow.enable.dsl = 2

// import modules
include { SplitGeneVariantPairs; ExtractVariants; ExtractGeneVariantPairs; ExtractLociBed; ExtractLociAll } from './modules/CollectResults'
include { AnnotateLoci } from './modules/CollectSignificantLoci'



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
--output         	  Path to outputfolder

""".stripIndent()

}

params.maf_table = 'NO_FILE'
params.bed = 'NO_FILE'
params.variants = 'NO_FILE'
params.gene_variant_pairs = 'NO_FILE'
params.inclusion_step_output = 'NO_FILE'
params.cols = '+p_value,+z_score'
params.p_threshold = 'NULL'
params.output

if (params.help){
    helpmessage()
    exit 0
}

//Default parameters
Channel.fromPath(params.input).collect().set { input_parquet_ch }
Channel.fromPath(params.genes).splitCsv(header: true).map { row -> "${row.ID}" } .set { genes_ch }
Channel.fromPath(params.variant_reference).collect().set { variant_reference_ch }
Channel.fromPath(params.gene_reference).collect().set { gene_reference_ch }

cohorts_ch = Channel.fromPath(params.mastertable)
    .ifEmpty { error "Cannot find master table from: ${params.mastertable}" }
    .splitCsv(header: true, sep: '\t', strip: true)
    .map{row -> [ row.cohort_new_name ]}
    .collect()

inclusion_dir_ch = file(params.inclusion_step_output)
bed_file_ch = file(params.bed)
variants_ch = file(params.variants)
gene_variant_pairs_ch = file(params.gene_variant_pairs)

Channel.fromPath(params.maf_table).collect().set { maf_table_ch }

gene_chunk_size=100
locus_chunk_size=100

extract_variants = params.variants != "NO_FILE"
extract_gene_variant_pairs = params.gene_variant_pairs != "NO_FILE"
extract_loci = params.bed != "NO_FILE"
annotate = params.maf_table != "NO_FILE"

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
summary['Permuted eQTLs']                           = params.permuted
summary['Reference data']                           = params.reference_data
summary['Genome reference']                         = params.genome_reference
summary['Variant reference']                        = params.variant_reference
summary['Gene reference']                           = params.gene_reference
summary['Gene list']                                = params.genes
summary['Extract loci']                             = extract_loci
summary['Extract variants']                         = extract_variants

log.info summary.collect { k,v -> "${k.padRight(21)}: $v" }.join("\n")
log.info "======================================================="

// Procedure:
// Generate breakpoints to get equal numbers of variants per bin
// Get uncorrelated variants per maf bin
// Using the uncorrelated variants, do accurate permutation p-value calculation

workflow {
    // Buffer genes
    genes_buffered_ch = genes_ch.collate(gene_chunk_size)

    if ( extract_gene_variant_pairs ) {
        SplitGeneVariantPairs(gene_variant_pairs_ch,gene_chunk_size)
        // Extract variants
        variants_extracted_ch = ExtractGeneVariantPairs(input_parquet_ch, variant_reference_ch, SplitGeneVariantPairs.out.flatten().view(), params.cols)
            .flatten()
            .map { file ->
                   def key = variants_ch.baseName
                   return tuple(key, file) }
            groupTuple()

        // Annotate loci
        variants_annotated_ch = AnnotateLoci(variants_extracted_ch, variant_reference_ch, gene_reference_ch, maf_table_ch, inclusion_dir_ch, cohorts_ch)

    }

    if ( extract_variants ) {
        // Extract variants
        variants_extracted_ch = ExtractVariants(input_parquet_ch, variant_reference_ch, genes_buffered_ch, variants_ch, params.cols)
            .flatten()
            .map { file ->
                   def key = variants_ch.baseName
                   return tuple(key, file) }
            groupTuple()

        // Annotate loci
        variants_annotated_ch = AnnotateLoci(variants_extracted_ch, variant_reference_ch, gene_reference_ch, maf_table_ch, inclusion_dir_ch, cohorts_ch)

    }

    if ( extract_loci ) {
        // Chunk loci
        bed_file_ch.splitText( by: locus_chunk_size )

        // Extract loci
        loci_extracted_ch = ExtractLociAll(input_parquet_ch, loci_ch, variant_reference_ch, genes_buffered_ch, params.cols)
            .flatten()
            .map { file ->
                   def key = file.name.toString().tokenize('.').get(1)
                   return tuple(key, file) }
            groupTuple()

        // Annotate loci
        loci_annotated_ch = AnnotateLoci(loci_extracted_ch, variant_reference_ch, gene_reference_ch, maf_table_ch, inclusion_dir_ch, cohorts_ch)

    }

    if ( extract_loci == false & extract_variants == false & extract_gene_variant_pairs == false) {
        // Extract all
        all_extracted_ch = ExtractVariants(input_parquet_ch, variant_reference_ch, genes_buffered_ch, variants_ch, params.cols, params.p_threshold)
            .flatten()
            .map { file ->
                   def key = file.name.toString().tokenize('.').get(1)
                   return tuple(key, file) }
            .groupTuple().view()

        // Annotate loci
        if (annotate) {
            variants_annotated_ch = AnnotateLoci(all_extracted_ch, variant_reference_ch, gene_reference_ch, maf_table_ch, inclusion_dir_ch, cohorts_ch)
        } else {
            all_extracted_ch
                .map{ row -> row[1] }
                .collectFile(name: 'extracted_merged.txt', skip: 1, keepHeader: true, storeDir: params.output)
        }

    }

}

workflow.onComplete {
    println ( workflow.success ? "Pipeline finished!" : "Something crashed...debug!" )
}
