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
        val loci

    output:
        path "ld.txt"

    script:
        """
        loci > "loci.txt"

        ld_calculator.py \
        --permuted ${permuted} \
        --genes ${genes} \
        --loci loci.txt \
        --out "ld.txt"
        """
}
