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
        variants_arg = (params.variants != '') ? "--variants-file ${variants}" : ""
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
            !{variants_arg} \
            --variant-reference !{variant_reference} \
            --output-prefix extracted \
            --cols '${cols}'

        rm -r tmp_eqtls
        '''
}


process ExtractLociBed {
    scratch true

    input:
        path input
        val loci
        path variantReference
        val genes
        val cols

    output:
        path "extracted*out.csv.gz"

    shell:
        gene_arg = genes.join(" ")
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
            --cols '!{cols}' \
            --bed-file "!{loci}" \
            --output-prefix extracted


        '''
}


process ExtractLociAll {
    scratch true

    input:
        path input
        val loci
        path variantReference
        val genes

    output:
        path "extracted*out.csv.gz"

    shell:
        gene_arg = genes.join(" ")
        phenotypes_formatted = genes.collect { "phenotype=$it" }.join("\n")
        '''
        mkdir tmp_eqtls
        echo "!{phenotypes_formatted}" > file_matches.txt

        while read gene; do
          cp -r "!{input}/${gene}" tmp_eqtls/
        done <file_matches.txt

        extract_parquet_results.py \
            --input-file tmp_eqtls \
            --genes !{gene_arg} \
            --variant-reference !{variantReference} \
            --cols '!{cols}' \
            --bed-file "!{loci}" \
            --output-prefix extracted

        '''
}