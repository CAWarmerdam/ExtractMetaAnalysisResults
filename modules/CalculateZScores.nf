#!/bin/bash nextflow


process CalculateZScores {
    scratch true

    input:
        path input
        path variant_reference
        val genes
        val variants

    output:
        path "z_scores.out.csv"

    shell:
        phenotypes_formatted = genes.collect { "phenotype=$it" }.join("\n")
        '''
        mkdir tmp_eqtls
        echo "!{phenotypes_formatted}" > file_matches.txt

        while read gene; do
          cp -r "!{input}/${gene}" tmp_eqtls/
        done <file_matches.txt

        extract_parquet_results.py \
            --input-file tmp_eqtls \
            --genes !{genes.join(' ')} \
            --variants-file !{variants.join(' ')} \
            --variant-reference !{variant_reference} \
            --output-prefix z_scores \
            --cols '+z_score'

        rm -r tmp_eqtls
        '''
}
