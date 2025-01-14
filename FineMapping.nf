#!/usr/bin/env nextflow

/*
 * enables modules
 */
nextflow.enable.dsl = 2

// import modules
include { AnnotateSignificantVariants; ExtractSignificantResults; DefineFineMappingLoci } from './modules/CollectSignificantLoci'
include { RunFineMappingOnCalculatedLd; UncorrelatedGenes; ExportResults } from './modules/RunFineMapping'
include { GetUncorrelatedVariants } from './modules/UncorrelatedVariants'

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

params.inclusion_step_output = 'NO_FILE'

if (params.help){
    helpmessage()
    exit 0
}


// Define list of chromosomes to analyse
chromosomes = params.chromosome ? [params.chromosome.toString()] : (1..22).collect { it.toString() }

//Default parameters
Channel.fromPath(params.significant).collect().set { significant_subset_ch }
Channel.fromPath(params.empirical).collect().set { empirical_parquet_ch }
Channel.fromPath(params.permuted).collect().set { permuted_parquet_ch }
Channel.fromPath(params.genes).splitCsv(header: ['gene']).map { row -> "${row.gene}" } .set { genes_ch }
Channel.fromPath(params.uncorrelated_genes).collect().set { uncorrelated_genes_ch }
Channel.fromPath(params.genome_reference).collect().set { genome_ref_ch }
Channel.fromPath(params.variant_reference).collect().set { variant_reference_ch }
Channel.fromPath(params.gene_reference).collect().set { gene_reference_ch }

gene_chunk_size=200
loci_per_job=1 // TODO: TMP, change back to 5

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
summary['Significant eQTLs']                        = params.significant
summary['Genome reference']                         = params.genome_reference
summary['Variant reference']                        = params.variant_reference
summary['Gene reference']                           = params.gene_reference
summary['Gene list']                                = params.genes
summary['Chromosomes']				    = chromosomes

log.info summary.collect { k,v -> "${k.padRight(21)}: $v" }.join("\n")
log.info "======================================================="

workflow LOCI {
    take:
        significant_subset_ch
        variant_reference_ch
        genome_ref_ch

    main:
        // Annotate significant variants with bp positions
        significant_results_ch = AnnotateSignificantVariants(significant_subset_ch, variant_reference_ch)
            .flatten()
            .map { file ->
             // Extract the chromosome number using regex
             def match = file.name =~ /chr(\w+)\.csv/
             def chromosome = match.size() == 1 ? match[0][1] : null
             // Return a tuple (chromosome number, file)
             [chromosome, file]
        }
        .filter { tuple -> tuple[0] != null } // Ensure chromosome is found
        .filter { tuple -> tuple[0] in chromosomes } // Filter on chromosome to analyse
        .view()

        // Output a channel of sorted bed files (by start pos).
        // Each bed file has the following columns: chr, start, end, name.
        // Where, each row is defined by a significant lead variant, with a 1Mb base pair window around the lead variant
        // The highest end pos must never be greater than the lowest start pos + 5Mb
        loci_ch = DefineFineMappingLoci(significant_results_ch, genome_ref_ch)
    emit:
        loci = loci_ch
}

workflow FINEMAPPING {
    take:
        empirical_parquet_ch
        permuted_parquet_ch
        variant_reference_ch
        uncorrelated_genes_ch
        loci_bed_ch
        loci_per_job

    main:
        loci_bed_collated_ch = loci_bed_ch.flatMap { chromosome, files -> files.collate(loci_per_job).collect { chunk -> [chromosome, chunk] } }

        finemapped_split_ch = RunFineMappingOnCalculatedLd(empirical_parquet_ch, permuted_parquet_ch, variant_reference_ch, uncorrelated_genes_ch, loci_bed_collated_ch).flatten().collect()

        // Write out results
        ExportResults(finemapped_split_ch)

}

workflow {
    // Buffer genes
    genes_buffered_ch = genes_ch.collate(gene_chunk_size)

    // Define loci to do finemapping for
    LOCI(
        significant_subset_ch,variant_reference_ch,genome_ref_ch)

    // Do finemapping
    FINEMAPPING(
        empirical_parquet_ch,permuted_parquet_ch,variant_reference_ch,
        uncorrelated_genes_ch.collect(), LOCI.out.loci, loci_per_job)
}

workflow.onComplete {
    println ( workflow.success ? "Pipeline finished!" : "Something crashed...debug!" )
}
