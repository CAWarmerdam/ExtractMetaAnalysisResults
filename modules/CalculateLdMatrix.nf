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
        cat !{uncorrelatedGenes} | sort | uniq > unique_genes_permuted.txt

        # Get set of genes to use in empirical analysis
        awk -F'\t' 'BEGIN {OFS = FS} {print $4}' !{bedFile} | sort | uniq > unique_genes_empirical.txt

        while read gene; do
          cp -r "!{empirical}/phenotype=${gene}" tmp_empirical/
        done <unique_genes_empirical.txt

        while read gene; do
          cp -r "!{permuted}/phenotype=${gene}" tmp_permuted/
        done <unique_genes_permuted.txt

        bedtools merge -i !{bedFile} -d 3000000 > ld_window.bed

        ld_calculator.py \
        --input-file "tmp_permuted" \
        --variant-reference !{variantReference} \
        --genes-file !{uncorrelatedGenes} \
        --bed-file ld_window.bed \
        --output-prefix "ld"

        cat finemapping_loci_cluster_2.bed | tr '\t'  ',' | while IFS=',' read -r chrom start end gene cluster; do
            echo "${chrom}\t${start}\t${end}\t${gene}\n" > "current_locus_as_bed_file.bed";
            echo "Extracting associations for ${chrom}:${start}-${end} and gene ${gene}";

        done
        '''
}
