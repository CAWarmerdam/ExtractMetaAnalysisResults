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
            --cols "p_value" \
            --output-file loci.txt
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
        annotate_loci.py \
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
        path bedFile
        path genomeRef

    output:
        path "merged.bed"

    script:
        // Define background bed file to take into account
        def bed = bedFile.name != 'NO_FILE' ? "$bedFile" : ''

        // Calculate flanks for genes, calculate flanks for snps, calculate union.
        """
        bedtools slop -i "${variantLoci}" -g "${genomeRef}" -b "${variantFlankSize}" > "variant_loci.flank.bed"
        bedtools slop -i "${geneLoci}" -g "${genomeRef}" -b "${geneFlankSize}" > "gene_loci.flank.bed"

        cat "variant_loci.flank.bed" "gene_loci.flank.bed" ${bed} > "total.flank.bed"

        # Get the union of the two bed files (including flanks)
        bedtools sort -i "total.flank.bed" > "total.flank.sorted.bed"
        bedtools merge -i "total.flank.sorted.bed" -d 0 -c 4 -o distinct > "merged.bed"
        """
}
