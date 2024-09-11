#!/usr/bin/env nextflow

/*
 * enables modules
 */
nextflow.enable.dsl = 2

// import modules


def helpmessage() {

log.info"""

HASE output analyzer v${workflow.manifest.version}
==============================================
Pipeline for parallelized extraction and filtering of the raw HASE results.
""".stripIndent()

}

params.maf_table = 'NO_FILE'
params.background_bed = 'NO_FILE'
params.inclusion_step_output = 'NO_FILE'
params.mastertable = 'NO_FILE'

if (params.help){
    helpmessage()
    exit 0
}

//Default parameters
Channel.fromPath(params.meta).collect().set { empirical_parquet_ch }
Channel.fromPath(params.cohort).collect().set { empirical_parquet_per_cohort_ch }
Channel.fromPath(params.genes).splitCsv(header: ['gene']).map { row -> "${row.gene}" } .set { genes_ch }
Channel.fromPath(params.genome_reference).collect().set { genome_ref_ch }
Channel.fromPath(params.variant_reference).collect().set { variant_reference_ch }
Channel.fromPath(params.gene_reference).collect().set { gene_reference_ch }

cohorts_ch = Channel.fromPath(params.mastertable)
    .ifEmpty { error "Cannot find master table from: ${params.mastertable}" }
    .splitCsv(header: true, sep: '\t', strip: true)
    .map{row -> [ row.cohort_new_name ]}
    .collect()

inclusion_step_output_ch = file(params.inclusion_step_output)
bed_file_ch = file(params.background_bed)

Channel.fromPath(params.maf_table).collect().set { maf_table_ch }

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

log.info summary.collect { k,v -> "${k.padRight(21)}: $v" }.join("\n")
log.info "======================================================="

workflow {
    // Buffer genes
    genes_buffered_ch = genes_ch.collate(gene_chunk_size)


}

workflow.onComplete {
    println ( workflow.success ? "Pipeline finished!" : "Something crashed...debug!" )
}
