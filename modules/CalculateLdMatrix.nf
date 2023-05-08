#!/bin/bash nextflow


process UncorrelatedGenes {

    publishDir "${params.output}/gene_correlations", mode: 'copy', overwrite: true

    input:
        path matrix
        val threshold

    output:
        path "uncorrelated_genes.txt", emit: genes
        path "gene_correlation_matrix.csv.gz", emit: correlations

    script:
        """
        export PYTHONUNBUFFERED="1"
        uncorrelated_genes.py ${matrix} "uncorrelated_genes.txt" -t ${threshold}
        """
}

process CalculateLdMatrix {

    publishDir "${params.output}/ld_matrices", mode: 'copy', overwrite: true

    input:
        val loci
        path results

    output:
        path "ld.*.csv.gz"

    shell:
        '''
        echo "${loci}" > "loci.txt"

        ld_calculator.py \
        --input-file ${results} \
        --loci loci.txt \
        --output-prefix "ld"
        '''
}
