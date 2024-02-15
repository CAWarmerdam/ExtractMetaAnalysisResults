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
        path "sign_variants.csv", optional:true

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
            --out-prefix annotated.

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
        path signVariants
        path genomeRef

    output:
        path "finemapping_loci_*.bed"

    shell:
        '''
        # Input is expected to be a table of significant results, with the following columns:
        # 7: Chromosome of the variant of the significant eQTL
        # 10: Basepair position of the variant of the significant eQTL
        # 2: ENSG identifier of the gene of the significant eQTL

        # First, get the relevant columns to make a bed file, and apply a splop to make a total
        # window of 3Mb
        awk -F'\t' 'BEGIN {OFS = FS} NR>1 {print $12,$9,$9,$2}' !{signVariants} \
        | bedtools slop -b 1500000 -g !{genomeRef} > loci_3Mb_window.bed

        # Second, extract the gene ENSG identifiers, and get a set of identifiers, removing duplicates
        awk -F'\t' 'BEGIN {OFS = FS} NR>1 {print $2}' !{signVariants} | sort | uniq > unique_genes.txt

        # Loop through the set of ENSG identifiers, merging the loci for each gene
        while read g; do
            grep "$g" loci_3Mb_window.bed | bedtools sort | bedtools merge -c 4 -o distinct >> loci_3Mb_merged_per_gene.bed
        done < unique_genes.txt

        # Sort the resulting bed file
        bedtools sort -i loci_3Mb_merged_per_gene.bed > loci_3Mb_merged_per_gene_sorted.bed

        # Assign your windows to clusters.
        assign_clusters.py loci_3Mb_merged_per_gene_sorted.bed 5000000 finemapping_loci
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
