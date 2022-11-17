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

""".stripIndent()

}

if (params.help){
    helpmessage()
    exit 0
}

//Default parameters
Channel.fromPath(params.inputfolder).set {parquet_path_ch}
Channel.fromPath(params.outputfile).set {output_file_ch}


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

process GenerateSqliteScript {
    tag {GenerateSqliteScript}

    input:
        path parquet_path from parquet_path_ch

    output:
        file("sqlite_script.sql") into sqlite_script_ch

    script:
        """
        python3 $baseDir/bin/generate_sqlite_script.py ${parquet_path}
        """
}

process GenerateSqliteDatabase {
    tag {GenerateSqliteDatabase}

    input:
        file sqlite_script from sqlite_script_ch
        file output_file from output_file_ch

    script:
        """
        sqlite3 ${output_file} '.read ${sqlite_script}'
        """
}

workflow.onComplete {
    println ( workflow.success ? "Pipeline finished!" : "Something crashed...debug!" )
}
