#!/bin/bash nextflow


process ExtractVariants {
    scratch true

    input:
        path input
        path variant_reference
        val genes
        path variants
        val cols

    output:
        path "extracted*.out.csv"

    shell:
        variants_arg = (variants.name != 'NO_FILE') ? "--variants-file ${variants}" : ""
        phenotypes_formatted = genes.collect { "phenotype=$it" }.join("\n")
        prefix = (genes.size() == 1) ? genes.collect { "extracted.$it" }.join("") : "extracted"
        '''
        mkdir tmp_eqtls
        echo "!{phenotypes_formatted}" > file_matches.txt

        while read gene; do
          cp -r "!{input}/${gene}" tmp_eqtls/
        done <file_matches.txt

        extract_parquet_results.py \
            --input-file tmp_eqtls \
            --genes !{genes.join(' ')} \
            !{variants_arg} \
            --variant-reference !{variant_reference} \
            --output-prefix !{prefix} \
            --cols '!{cols}'

        rm -r tmp_eqtls
        '''
}


process ExtractLociBed {
    scratch true

    input:
        path input
        path locus
        path variantReference
        val genes
        val cols

    output:
        path "extracted*out.csv"

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
            --variant-reference !{variantReference} \
            --genes !{genes.join(' ')} \
            --cols '!{cols}' \
            --bed-file !{locus} \
            --output-prefix extracted


        '''
}


process ExtractLociAll {
    scratch true

    input:
        path input
        path locus
        path variantReference
        val genes
        val cols

    output:
        path "extracted*out.csv"

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
            --variant-reference !{variantReference} \
            --genes !{genes.join(' ')} \
            --cols '!{cols}' \
            --bed-file !{locus} \
            --output-prefix extracted
        '''
}