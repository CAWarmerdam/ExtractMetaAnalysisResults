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


process FindLdPartners {
    //scratch true // Needs to be set to true!
    scratch '$TMPDIR'
    publishDir "${params.output}/ld_partners", mode: 'copy', overwrite: true

    input:
        path variant_path
        path ld_panel, stageAs: 'ld_panel'
        path variantReference

    output:
        path "*_ld_partners.txt.gz", emit: output

    shell:
        '''
        find_ld_partners.R \
          --independent-variants !{variant_path} \
          --ld !{ld_panel} \
          --variant-reference !{variantReference} \
          --ld-type pcs \
        '''
}


//Default parameters
Channel.fromPath(params.independent_variants).collect().splitText( by: 1000 ).set { variant_file_ch }
Channel.fromPath(params.ld).collect().set { permuted_parquet_ch }
Channel.fromPath(params.variant_reference).collect().set { variant_reference_ch }

workflow {
    FindLdPartners(variant_file_ch, permuted_parquet_ch, variant_reference_ch)
}

workflow.onComplete {
    println ( workflow.success ? "Pipeline finished!" : "Something crashed...debug!" )
}
