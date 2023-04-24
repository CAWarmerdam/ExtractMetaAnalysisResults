#!/bin/bash nextflow


process CalculateZScores {

    input:
        path input
        val genes
        val variants

    output:
        path "z_scores.txt"

    script:
        """
        python3 $baseDir/bin/extract_parquet_results.py \
            --input-file ${input} \
            --genes ${genes.join(' ')} \
            --variants ${variants.join(' ')} \
            --output-file z_scores.txt \
            --cols 'z_scores' \
        """
}
