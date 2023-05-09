#!/bin/bash nextflow


process SampleOverlapMatrix {

    publishDir "${params.output}/sample_overlap", mode: 'copy', overwrite: true

    input:
        path inclusionDir

    output:
        path "sample_overlap_matrix.csv"

    script:
        """
        sample_overlap_table.py ${inclusionOutput}/filter_logs_full.log ${inclusionOutput}/filter_logs_full.log
        """
}


process RunHyprColoc {

    input:
        tuple val(locus_string), path(empirical), path(permuted)
        path geneCorrelations
        path sampleOverlap
        val posteriorThreshold
        val csThreshold
        val outputCsPip

    output:
        path "cis_trans_coloc.cs.txt"
        path "cis_trans_coloc.variant_pips.txt"

    script:
    """
    Rscript --vanilla $baseDir/bin/RunHyprColoc.R \
    --eqtls ${eQtlResults} \
    --output-cs-pip ${outputCsPip} \
    --posterior-threshold ${posteriorThreshold} \
    --cs-threshold ${csThreshold} \
    --output-prefix "cis_trans_coloc"
    """
    if (outputCsPip == "FALSE")
      """
      touch cis_trans_coloc.variant_pips.txt
      """
}