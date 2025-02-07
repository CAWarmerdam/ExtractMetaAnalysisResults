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

process AnnotateSignificantVariants {

    input:
        path signSubset
        path variantReference

    output:
        path "sign_variants.chr*.csv"

    shell:
        '''
        annotate_significant_variants.R --sign-subset !{signSubset} --variant-reference !{variantReference} --output-prefix "sign_variants"
        '''
}

process DefineFineMappingLoci {
    input:
        tuple val(chromosome), path(signVariants)
        path genomeRef
        val windowMb

    output:
        tuple val(chromosome), path("finemapping_loci_*.bed")

    shell:
        '''
        # Input is expected to be a table of significant results, with the following columns:
        # 7: Chromosome of the variant of the significant eQTL
        # 10: Basepair position of the variant of the significant eQTL
        # 2: ENSG identifier of the gene of the significant eQTL

        # First, get the relevant columns to make a bed file, and apply a splop to make a total
        # window of 3Mb
        awk -F'\t' 'BEGIN {OFS = FS} NR>1 {print $13,$12-1,$12,$1}' !{signVariants} \
        | bedtools slop -b !{windowMb/2*1000000} -g !{genomeRef} > loci_3Mb_window.bed

        # Second, extract the gene ENSG identifiers, and get a set of identifiers, removing duplicates
        awk -F'\t' 'BEGIN {OFS = FS} NR>1 {print $1}' !{signVariants} | sort | uniq > unique_genes.txt

        # Loop through the set of ENSG identifiers, merging the loci for each gene
        while read g; do
            grep "$g" loci_3Mb_window.bed | bedtools sort | bedtools merge -c 4 -o distinct >> loci_3Mb_merged_per_gene.bed
        done < unique_genes.txt

        # Potentially remove the HLA region
        echo '6\t28510120\t33480577\tHLA\n' > hla_range.bed
        bedtools intersect -a loci_3Mb_merged_per_gene.bed -b hla_range.bed -v > loci_3Mb_merged_per_gene_filtered.bed

        # Sort the resulting bed file
        bedtools sort -i loci_3Mb_merged_per_gene_filtered.bed > loci_3Mb_merged_per_gene_filtered_sorted.bed

        # Assign your windows to clusters.
        assign_clusters.py loci_3Mb_merged_per_gene_filtered_sorted.bed 5000000 finemapping_loci
        '''
}

process AnnotateLoci {
    publishDir "${params.output}/loci_empirical_annotated", mode: 'copy', overwrite: true

    input:
        path "locus_*.csv"
        path variantReference
        path geneReference
        path mafTable
        path inclusionDir
        val cohorts

    output:
        tuple val(locus_string), path("annotated.csv.gz")

    script:
        """
        annotate_loci.py \
            --input-file locus_*.csv \
            --cohorts ${cohorts.join(' ')} \
            --variant-reference ${variantReference} \
            --gene-gff ${geneReference} \
            --maf-table ${mafTable} \
            --inclusion-path ${inclusionDir} \
            --out-prefix annotated
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
