#!/bin/bash nextflow


process SplitGeneVariantPairs {
    executor 'local'

    input:
        path gene_variant_file
        path genes
        val n_genes_per_chunk

    output:
        path "output_chunk_*.csv"

    shell:
        '''
        gunzip !{gene_variant_file} -c > gene_variant_file.txt
        head -n 1 gene_variant_file.txt > gene_variant_file_pruned.txt
        grep -Ff !{genes} gene_variant_file.txt >> gene_variant_file_pruned.txt
        split_variant_gene_pairs.py -P gene_variant_file_pruned.txt -k !{n_genes_per_chunk}
        '''
}


process ExtractVariants {
    scratch true
    publishDir "${params.output}", mode: 'copy', overwrite: true, enabled: true, saveAs: { fn -> "${genes[0]}.out.csv" }

    input:
        path input
        path variant_reference
        val genes
        path variants
        val cols
        val p_threshold
        val mode

    output:
        path "extracted*.out.csv"

    shell:
        variants_arg = (variants.name != 'NO_FILE') ? "--variants-file ${variants}" : ""
        p_threshold_arg = (p_threshold != 'NULL') ? "--p-thresh ${p_threshold}" : ""
        clump_arg = (mode == 'clump') ? "--dist-clump" : ""
        phenotypes_formatted = genes.collect { "phenotype=$it" }.join("\n")
        prefix = (genes.size() == 1) ? genes.collect { "extracted.$it" }.join("") : "extracted.${genes[0]}"
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
            !{p_threshold_arg} \
            !{clump_arg} \
            --variant-reference !{variant_reference} \
            --output-prefix !{prefix} \
            --cols '!{cols}'

        rm -r tmp_eqtls
        '''
}


process ExtractCorrectedTransQtls {
    scratch true

    input:
        path input
        path cisExplainedVariance
        path variantReference
        path geneReference
        val genes

    output:
        path "extracted*.txt"

    shell:
        phenotypes_formatted = genes.collect { "phenotype=$it" }.join("\n")
        genes_formatted = genes.join("\n")
        '''
        mkdir tmp_eqtls
        echo "!{phenotypes_formatted}" > file_matches.txt
        echo "!{genes_formatted}" > genes.txt

        while read gene; do
          cp -rL "!{input}/${gene}" tmp_eqtls/
        done <file_matches.txt

        correct_trans_for_cis.R \
            --input-path tmp_eqtls \
            --cis-explained-variance !{cisExplainedVariance} \
            --variant-reference !{variantReference} \
            --gene-reference !{geneReference} \
            --genes genes.txt \
            --output-prefix extracted

        rm -r tmp_eqtls
        '''
}


process ExtractGeneVariantPairs {
    scratch true

    input:
        path input
        path variant_reference
        path gene_variant_pairs
        val cols

    output:
	    path "*.extracted.out.csv"

    shell:
        out_prefix = gene_variant_pairs.getSimpleName()
        '''
        mkdir tmp_eqtls
        awk '{FS="\t"; OFS="\t"} {print "phenotype="$2}' !{gene_variant_pairs} | sort | uniq > file_matches.txt

        while read gene; do
          [[ -e "!{input}/${gene}" ]] && cp -r "!{input}/${gene}" tmp_eqtls/
        done <file_matches.txt

        extract_parquet_results.py \
            --input-file tmp_eqtls \
            -P !{gene_variant_pairs} \
            --variant-reference !{variant_reference} \
            --output-prefix "!{out_prefix}.extracted" \
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
