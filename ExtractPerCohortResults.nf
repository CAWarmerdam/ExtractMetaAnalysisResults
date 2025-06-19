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
    input:
        master_table
        input_parquet
        variant_reference
        genes
        maf_table
        finemapped_variants

    output:
        "extracted_finemapped_variants_per_cohort_*.txt.gz"

    shell:
        '''
        extract_per_cohort_results.R \
        --input-parquet !{input_parquet} \
        --master-table !{master_table} \
        --gene_list !{genes.join(" ")} \
        --maf-table !{input_maf} \
        --variant-reference !{variant_reference} \
        --finemapped-variants !{finemapped_variants}
        '''
}

workflow {
    // Buffer genes
    genes_buffered_ch = genes_ch.collate(gene_chunk_size)

    per_cohort_results_ch = ExtractPerCohortResults(master_table_ch, input_parquet_ch, variant_reference_ch, genes_buffered_ch, maf_table_ch, finemapped_variants_ch)

    per_cohort_results_ch.collect()
}

workflow.onComplete {
    println ( workflow.success ? "Pipeline finished!" : "Something crashed...debug!" )
}
