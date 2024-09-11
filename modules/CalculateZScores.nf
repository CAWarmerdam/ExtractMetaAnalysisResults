#!/bin/bash nextflow


process CalculateZScoreMatrix {
    scratch true

    input:
        path input
        path variant_reference
        val genes
        val variants
        val n_threshold
        val cohorts

    output:
        path "z_scores.out.z_score.csv", emit: z_scores
        path "z_scores.out.sample_size.csv", emit: sample_size

    shell:
        phenotypes_formatted = genes.collect { "phenotype=$it" }.join("\n")
        cohort_arg = (cohorts.size() == 1) ? "--cohort ${cohorts[0]}" : ""
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
            --n-threshold !{n_threshold} \
            --as-matrix \
            --cols 'sample_size,+z_score' !{cohort_arg}

        rm -r tmp_eqtls
        '''
}

process ConcatMatrix {
    publishDir "${params.output}/exported_matrix", mode: 'copy', overwrite: true

    input:
        path matrix, stageAs: "z_scores.out.*.csv"
        val name

    output:
        path "${name}"

    shell:
        '''
        concat_matrix.py \
            --input !{matrix.join(' ')} \
            --output !{name}
        '''
}
