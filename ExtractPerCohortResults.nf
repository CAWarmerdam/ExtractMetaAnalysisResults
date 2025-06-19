#!/usr/bin/env nextflow

/*
 * enables modules
 */
nextflow.enable.dsl = 2

def helpmessage() {

log.info"""

HASE output analyzer v${workflow.manifest.version}
==============================================
Pipeline for parallelized extraction and filtering of the per cohort raw HASE results.

This pipeline is used to extract subsets of results from the HASE results (numerous large .parquet files).

""".stripIndent()

}

process ExtractPerCohortResults {
    publishDir "${params.output}/extracted_summary_statistics_per_cohort", mode: 'copy', overwrite: true
    memory '8 GB'
    executor 'slurm'
    time '1h'

    input:
        path input_parquet
        path variant_reference
        path maf_table
        path finemapped_variants
        val genes

    output:
        "extracted_finemapped_variants_per_cohort_*.txt.gz"

    shell:
        '''
        extract_per_cohort_results.R \
        --input-parquet !{input_parquet} \
        --gene_list !{genes.join(" ")} \
        --maf-table !{input_maf} \
        --variant-reference !{variant_reference} \
        --finemapped-variants !{finemapped_variants}
        '''
}

//Default parameters
Channel.fromPath(params.input).collect().set { input_parquet_ch }
Channel.fromPath(params.genes).splitCsv(header: true).view().map { row -> "${row.ID}" } .set { genes_ch }
Channel.fromPath(params.variant_reference).collect().set { variant_reference_ch }
Channel.fromPath(params.maf_table).collect().set { maf_table_ch }
Channel.fromPath(params.eqtls).collect().set { finemapped_variants_ch }

gene_chunk_size=100

workflow {
    // Buffer genes
    genes_buffered_ch = genes_ch.collate(gene_chunk_size)

    ExtractPerCohortResults(input_parquet_ch, variant_reference_ch, maf_table_ch, finemapped_variants_ch, genes_buffered_ch)
}

workflow.onComplete {
    println ( workflow.success ? "Pipeline finished!" : "Something crashed...debug!" )
}
