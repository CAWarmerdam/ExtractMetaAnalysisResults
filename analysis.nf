#!/usr/bin/env nextflow

def helpmessage() {

log.info"""

HASE output analyzer v${workflow.manifest.version}
==============================================
Pipeline for parallelized extraction and filtering of the raw HASE results.

This pipeline is used to extract subsets of results from the HASE results (numerous large .parquet files).

Usage:

nextflow run ExtractHaseResults.nf \
--inputfolder '/inputfolder/' \
--outputfile '/outputfile/' \

Mandatory arguments:
--inputfolder           Path to the folder with HASE result .parquet files.
--outputfile            Path to where the database should be written
--outputfolder		Path to output

""".stripIndent()

}

if (params.help){
    helpmessage()
    exit 0
}

//Default parameters
Channel.fromPath(params.inputfolder).collect().set { parquet_path_ch }


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
summary['Input directory']                          = params.inputfolder
summary['Output file']                              = params.outputfile

log.info summary.collect { k,v -> "${k.padRight(21)}: $v" }.join("\n")
log.info "======================================================="

process Analysis {
    tag {Analysis}

    publishDir "${params.outputfolder}", mode: 'copy', overwrite: true

    input:
        path parquet_dataset from parquet_path_ch.collect()
        path variant_set from Channel.fromPath(params.variants)

    output:
        file("ExtractedPhenotype_*.txt.gz")

    script:
        """
        Rscript $baseDir/bin/db_query_template.R --dataset ${parquet_dataset} --variants ${variant_set}
        """

}

workflow.onComplete {
    println ( workflow.success ? "Pipeline finished!" : "Something crashed...debug!" )
}
