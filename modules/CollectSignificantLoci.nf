#!/bin/bash nextflow

process ExtractSignificantResults {
    scratch true

    input:
        path input
        path variantReference
        path geneReference
        path inclusionDir
        val genes
        val p_value
        val cohorts

    output:
        path "lead_variants.csv"

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
            --p-thresh !{p_value} \
            --cols "+p_value" \
            --output-prefix extracted

        annotate_lead_variants.py \
            --input-file extracted.out.csv \
            --cohorts !{cohorts.join(' ')} \
            --inclusion-path !{inclusionDir} \
            --variant-reference !{variantReference} \
            --gene-ref !{geneReference} \
            --out-prefix annotated.!{genes.join("_")}

        rm -r tmp_eqtls
        '''
}

process AnnotateResults {
    input:
        path significantResults
        path variantReference
        path geneReference

    output:
        path "loci.variants.bed"

    script:
        """
        eqtls_to_bed.py \
            --input-file ${significantResults} \
            --variant-reference ${variantReference} \
            --gene-ggf ${geneReference} \
            --out-prefix "loci"
        """
}

process DefineFineMappingLoci {
    input:
        path leadVariants
        path genomeRef

    output:
        path "finemapping_loci_*.bed"

    shell:
        '''
        awk '{ print $0,$1,$2,$3 }' !{leadVariants} | bedtools slop -b 1000000 -g !{genomeRef} > loci_1Mb_window.bed

        awk -F'\t' 'BEGIN {OFS = FS} NR>1 {print $7,$10,$10,$2}' !{leadVariants} \
        | bedtools slop -b 1000000 -g !{genomeRef} \
        | bedtools sort > loci_1Mb_window.bed

        # This is suboptimal since there will be ld chunks of over humongous size, use Dans method
        assign_clusters.py  loci_1Mb_window.bed 1000000 finemapping_loci
        '''
}

process AnnotateLoci {
    publishDir "${params.output}/loci_empirical_annotated", mode: 'copy', overwrite: true

    input:
        tuple val(locus_string), path(files, stageAs: "locus_*.csv")
        path variantReference
        path geneReference
        path mafTable
        path inclusionDir
        val cohorts

    output:
        tuple val(locus_string), path("annotated.${locus_string}.csv.gz")

    script:
        """
        head -n 1 ${files[0]} > concatenated.${locus_string}.csv
        tail -n +2 ${files.join(' ')} >> concatenated.${locus_string}.csv

        annotate_loci.py \
            --input-file concatenated.${locus_string}.csv \
            --cohorts ${cohorts.join(' ')} \
            --variant-reference ${variantReference} \
            --gene-gff ${geneReference} \
            --maf-table ${mafTable} \
            --inclusion-path ${inclusionDir} \
            --out-prefix annotated.${locus_string}
        """
}

process ConcatLoci {
    publishDir "${params.output}/loci_permuted", mode: 'copy', overwrite: true

    input:
        tuple val(locus_string), path(files, stageAs: "locus_*.csv")

    output:
        path "concatenated.${locus_string}.csv.gz"

    script:
        """
        head -n 1 ${files[0]} > concatenated.${locus_string}.csv
        tail -n +2 ${files.join(' ')} >> concatenated.${locus_string}.csv
        gzip -f concatenated.${locus_string}.csv
        """
}
