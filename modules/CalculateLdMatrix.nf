#!/bin/bash nextflow


process UncorrelatedGenes {
    input:
        path matrix
        val threshold

    output:
        path "uncorrelated_genes.txt"

    script:
        """
        uncorrelated_genes.py ${matrix} "uncorrelated_genes.txt" -t 0.1
        """
}

process CalculateLdMatrix {
    input:
        path permuted
        path genes
        path variantReference
        val loci

    output:
        path "ld.*.csv"

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
