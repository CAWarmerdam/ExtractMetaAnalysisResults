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
Channel.fromPath(params.inputfolder).into { parquet_path_ch; parquet_path_ch2 }


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
        echo ".load $baseDir/bin/libparquet" > sqlite_script.sql
        python3 $baseDir/bin/generate_sqlite_script.py ${parquet_path}/ >> sqlite_script.sql
        """
}

process GenerateSqliteDatabase {
    tag {GenerateSqliteDatabase}

    publishDir "${params.outputfolder}", mode: 'copyNoFollow', overwrite: true

    input:
        path parquet_path from parquet_path_ch2
        file sqlite_script from sqlite_script_ch
        val output_file from params.outputfile

    output:
        file("${output_file}")
        path("${parquet_path}")

    script:
        """
        sqlite3 ${output_file} '.read ${sqlite_script}'
        """
}

workflow.onComplete {
    println ( workflow.success ? "Pipeline finished!" : "Something crashed...debug!" )
}
