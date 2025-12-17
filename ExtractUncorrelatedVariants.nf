#!/usr/bin/env nextflow

/*
 * enables modules
 */
nextflow.enable.dsl = 2

// import modules
include { GetUncorrelatedVariants } from './modules/UncorrelatedVariants'


def helpmessage() {

log.info"""

HASE output analyzer v${workflow.manifest.version}
==============================================
Nextflow pipeline that selects uncorrelated variants
""".stripIndent()

}

if (params.help){
    helpmessage()
    exit 0
}

//Default parameters
Channel.fromPath(params.reference_data).set { reference_bcf_files_ch }

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

// Procedure:
// Generate breakpoints to get equal numbers of variants per bin
// Get uncorrelated variants per maf bin
// Using the uncorrelated variants, do accurate permutation p-value calculation

workflow {
    uncorrelated_variants_ch = GetUncorrelatedVariants(reference_bcf_files_ch)
            .collectFile(name: 'merged.prune.in', newLine: true, storeDir: params.outdir).collect()
}

workflow.onComplete {
    println ( workflow.success ? "Pipeline finished!" : "Something crashed...debug!" )
}
