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
Pipeline for extracting LD partners from a set of variants
Usage:
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
        path "*_ld_partners.txt", emit: output

    shell:
        '''
        find_ld_partners.R \
          --independent-variants !{variant_path} \
          --ld !{ld_panel} \
          --variant-reference !{variantReference} \
          --ld-type pcs \
          --r2-threshold 0.9 \
        '''
}


//Default parameters
Channel.fromPath(params.independent_variants).collect().splitText( by: 200, file:true ).set { variant_file_ch }
Channel.fromPath(params.ld).collect().set { permuted_parquet_ch }
Channel.fromPath(params.variant_reference).collect().set { variant_reference_ch }

workflow {
    FindLdPartners(variant_file_ch, permuted_parquet_ch, variant_reference_ch).collectFile(name: "ld_partners.txt", keepHeader:true, skip: 1, storeDir: params.output)
}

workflow.onComplete {
    println ( workflow.success ? "Pipeline finished!" : "Something crashed...debug!" )
}
