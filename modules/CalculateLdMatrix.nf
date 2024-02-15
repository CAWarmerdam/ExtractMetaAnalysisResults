#!/bin/bash nextflow


process UncorrelatedGenes {
    cache 'lenient'
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
        uncorrelated_genes.py --zscores-file ${matrix} --output-file "uncorrelated_genes.txt" -t ${threshold}
        """
}

process CalculateLdMatrix {
    scratch false
    publishDir "${params.output}/ld_matrices", mode: 'copy', overwrite: true

    input:
        path empirical, stageAs: 'empirical'
        path permuted, stageAs: 'permuted'
        path variantReference
        path uncorrelatedGenes
        path bedFile

    output:
        path "ld.*.csv.gz"

    shell:
        '''
        mkdir tmp_empirical
        mkdir tmp_permuted

        # Get set of unique genes to use in permuted analysis
        cat !{uncorrelatedGenes} | sort | uniq > unique_genes_empirical.txt

        # Get set of genes to use in empirical analysis
        awk -F'\t' 'BEGIN {OFS = FS} NR>1 {print $4}' !{bedFile} | sort | uniq > unique_genes_permuted.txt

        while read gene; do
          cp -r "!{empirical}/phenotype=${gene}" tmp_eqtls/
        done <unique_genes_empirical.txt

        while read gene; do
          cp -r "!{permuted}/phenotype=${gene}" tmp_eqtls/
        done <unique_genes_permuted.txt

        bedtools merge -i !{bedFile} -d 3000000 > ld_window.bed

        ld_calculator.py \
        --input-file "tmp_permuted" \
        --variant-reference !{variantReference} \
        --bed-file ld_window.bed \
        --output-prefix "ld"
        '''
}
