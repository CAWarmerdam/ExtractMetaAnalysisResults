#!/bin/bash nextflow

process ExtractSignificantResults {
    input:
        path eqtls
        val genes
        val p_value

    output:
        path "loci.txt"

    script:
        gene_arg = genes.join(" ")

        """
        extract_parquet_results.py \
            --input-file ${eqtls} \
            --genes ${gene_arg} \
            --p-thresh ${p_value} \
            --output-file z_scores.txt
        """
}

process AnnotateLoci {
    input:
        path significantResults
        path variantReference
        path geneReference

    output:
        path "loci.variants.bed", emit: variant_loci
        path "loci.genes.bed", emit: gene_loci

    script:
        """
        python3 $baseDir/bin/annotate_loci.py \
            --input-file ${significantResults} \
            --variant-reference ${variantReference} \
            --gene-ggf ${geneReference} \
            --out-prefix "loci"
        """
}

process IntersectLoci {
    input:
        path variantLoci
        val variantFlankSize
        path geneLoci
        val geneFlankSize
        path genomeRef

    output:
        path "union.bed"

    script:
        // Calculate flanks for genes, calculate flanks for snps, calculate union.
        """
        bedtools flank -i "${variantLoci}" -g "${genomeRef}" -b "${variantFlankSize}" > "variant_loci.flank.bed"
        bedtools flank -i "${geneLoci}" -g "${genomeRef}" -b "${geneFlankSize}" > "gene_loci.flank.bed"

        # Get the union of the two bed files (including flanks)
        bedtools unionbedg -i "variant_loci.flank.bed" "gene_loci.flank.bed" -d 0 > union.bed
        """
}
