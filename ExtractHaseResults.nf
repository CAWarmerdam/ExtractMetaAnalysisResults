#!/usr/bin/env nextflow

def helpmessage() {

log.info"""

HASE output analyzer v${workflow.manifest.version}
Pipeline for parallelized extraction and filtering of the raw HASE results.

This pipeline is used to extract subsets of results from the HASE results (numerous large .parquet files).

Usage:

nextflow run ExtractHaseResults.nf \
--inputfolder '/inputfolder/' \
--outputfile '/outputfile/' \
--NumberOfCovariates [number of covariates] \
--PhenFilter None \
--SnpFilter None \
--PThresh None \
--SnpRef '/SnpReferenceFromMapper/'
--Chunks [nr of chunks]

Mandatory arguments:
--inputfolder           Path to the folder with HASE result .parquet files.
--outputfile            Path to output file where results are written as tab-delimited text file.
--NumberOfCovariates    Number of covariates in the model. This is used for calculating the degrees of freedom and subsequently P-value.
--SnpRef                Reference file for variant mapper, as specified in HASE documentation. Has to be gzipped and space-delimited. This is used to include information about alleles.

Optional arguments:
--PThresh               P-value threshold for filtering the results. Defaults None.
--PhenFilter            Filter to include phenotypes into the association results. One of two: first, file with the list of phenotypes to include. In that case, no header expected. Second, individual phenotype IDs specified and separated by space. Defaults to None.
--SnpFilter             Filter to include SNPs/variants into the association results.  One of two: first, file with the list of SNPs/variants to include. In that case, no header expected. Second, individual SNP IDs specified and separated by space. Defaults to None.
--Chunks                Number of chunks which was used for original meta-analysis. Defaults to 100.
""".stripIndent()

}

if (params.help){
    helpmessage()
    exit 0
}

//Default parameters
params.PhenFilter = 'None'
params.SnpFilter = 'None'
params.Chunks = 100
params.PThresh = 1

InpFile = Channel.fromPath(params.inputfolder)
OutpFile = Channel.fromPath(params.outputfile)
NumberOfCovariates = Channel.value(params.NumberOfCovariates)
PhenFilter = Channel.fromPath(params.PhenFilter)
SnpFilter = Channel.fromPath(params.SnpFilter)
SnpRef = Channel.fromPath(params.SnpRef)


log.info """=======================================================
HASE output analyzer v${workflow.manifest.version}"
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
summary['SNP reference file']                       = params.SnpRef
summary['Number of covariates']                     = params.NumberOfCovariates
summary['P threshold']                              = params.PThresh
summary['Phenotype filter']                         = params.PhenFilter
summary['SNP filter']                               = params.SnpFilter
summary['Chunks']                                   = params.Chunks

log.info summary.collect { k,v -> "${k.padRight(21)}: $v" }.join("\n")
log.info "======================================================="

process ExtractResults {

    tag {"${Chunk}"}

    input:
        file InpFile from InpFile
        val NumberOfCovariates from NumberOfCovariates
        file PhenFilter from PhenFilter
        file SnpFilter from SnpFilter
        file SnpRef from SnpRef
        each Chunk from 1..params.Chunks
        val PThresh from params.PThresh

    output:
        file 'OutputChunk*.txt' into ProcessedOutput

    """
    nr_of_files=\$(ls -hlt ${InpFile}/node_${Chunk}_*.parquet 2>/dev/null | wc -l || true)

    echo "Analyzing \${nr_of_files} files"

    if [ "\$nr_of_files" -gt 0 ]
    then
        python2 $baseDir/bin/HaseOutputParquetAnalyzer.py \
        -i "${InpFile}/node_${Chunk}_*.parquet" \
        -o OutputChunk${Chunk}.txt \
        -nc ${NumberOfCovariates} \
        -sref ${SnpRef} \
        -phen ${PhenFilter} \
        -snps ${SnpFilter} \
        -p ${PThresh}
    else
        echo -e "phenotype\\tSNP\\tref_all\\talt_all\\tbeta\\tse\\tP\\ti_squared\\tN\\n" > OutputChunk${Chunk}.txt
    fi
    """

}

ProcessedOutput.collectFile(name: params.outputfile, keepHeader: true, sort: true)

workflow.onComplete {
    println ( workflow.success ? "Pipeline finished!" : "Something crashed...debug!" )
}
