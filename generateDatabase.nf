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

process MergeParquetFiles {
    tag {MergeParquetFiles}

    input:
        path parquet_path from parquet_path_ch
        val node_number from Channel.from( 1,2,3 )

    output:
        file("node_{node_number}.parquet") into parquet_combined_ch

    script:
        """
        python3 $baseDir/bin/coalesce_parquet.py \
        ${parquet_path} 'node_${node_number}_*.parquet' "node_{node_number}.parquet"
        """
}

process CollectParquetFiles {
    tag {CollectParquetFiles}

    input:
        file(parquet_files) from parquet_combined_ch.collect()

    output:
        path("data") into combined_parquet_path_ch_pre

    script:
        """
        mv "*.parquet" "data"
        """
}

combined_parquet_path_ch_pre.collect().set { combined_parquet_path_ch }

process GenerateSqliteScript {
    tag {GenerateSqliteScript}

    input:
        path parquet_path from combined_parquet_path_ch

    output:
        file("sqlite_script.sql") into sqlite_script_ch

    script:
        """
        echo ".load /tools/libparquet" > sqlite_script.sql
        python3 $baseDir/bin/generate_sqlite_script.py --path ${parquet_path} --parquet 'node_*.parquet' >> sqlite_script.sql
        """
}

process GenerateSqliteDatabase {
    tag {GenerateSqliteDatabase}

    publishDir "${params.dbfolder}", mode: 'copyNoFollow', overwrite: true

    input:
        path parquet_path from combined_parquet_path_ch
        file sqlite_script from sqlite_script_ch
        val output_file from params.dbfile

    output:
        file("${output_file}") into db_file_ch
        path("${parquet_path}")

    script:
        """
        sqlite3 ${output_file} '.read ${sqlite_script}'
        """
}

process Analysis {
    tag {Analysis}

    publishDir "{params.outputfolder}", mode: 'copy', overwrite: true

    input:
        path parquet_path from combined_parquet_path_ch
        file db_file from db_file_ch

    output:
        file("output.txt")

    script:
        """
        Rscript $baseDir/bin/db_query_template.R --database ${db_file}
        """

}

workflow.onComplete {
    println ( workflow.success ? "Pipeline finished!" : "Something crashed...debug!" )
}
