#!/bin/bash nextflow


process UncorrelatedGenes {

    publishDir "${params.output}", mode: 'copy', overwrite: true

    input:
        path matrix
        val threshold

    output:
        path "uncorrelated_genes.txt", emit: genes
        path "gene_correlation_matrix.csv.gz", emit: correlations

    script:
        """
        uncorrelated_genes.py ${matrix} "uncorrelated_genes.txt" -t 0.1
        """
}

process CalculateLdMatrix {

    publishDir "${params.output}", mode: 'copy', overwrite: true

    input:
        path permuted
        path genes
        path variantReference
        val loci

    output:
        path "ld.*.csv.gz"

    script:
        """
        echo "${loci}" > "loci.txt"

        ld_calculator.py \
        --input-file ${permuted} \
        --genes ${genes} \
        --loci loci.txt \
        --variant-reference ${variantReference} \
        --output-prefix "ld"
        """
}
