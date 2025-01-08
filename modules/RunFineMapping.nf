#!/bin/bash nextflow


process UncorrelatedGenes {
    publishDir "${params.output}/gene_correlations", mode: 'copy', overwrite: true

    input:
        path matrix_z
        path matrix_n
        val n_threshold
        val r_threshold

    output:
        path "uncorrelated_genes.txt", emit: genes
        path "gene_correlation_matrix.csv.gz", emit: correlations

    script:
        """
        export PYTHONUNBUFFERED="1"
        uncorrelated_genes.py --input-prefix 'mat' --output-file "uncorrelated_genes.txt" -n ${n_threshold} -t ${r_threshold}
        """
}

process RunFineMappingOnCalculatedLd {
    scratch true // Needs to be set to true!
    publishDir "${params.output}/finemapped", mode: 'copy', overwrite: true

    input:
        path empirical, stageAs: 'empirical'
        path permuted, stageAs: 'permuted'
        path variantReference
        path uncorrelatedGenes
        path bedFile

    output:
        path "finemapped.*.tsv"

    shell:
        '''
        mkdir tmp_empirical
        mkdir tmp_permuted

        # Get set of unique genes to use in permuted analysis
        cat !{uncorrelatedGenes} | sort | uniq > unique_genes_permuted.txt

        # Get set of genes to use in empirical analysis
        cat !{bedFile.join(" ")} > "all_loci.bed"
        awk -F'\t' 'BEGIN {OFS = FS} {print $4}' "all_loci.bed" | sort | uniq > unique_genes_empirical.txt

        # Need to add -L to cp command when running on a compute nodes scratch space
        while read gene; do
          cp -rL "!{empirical}/phenotype=${gene}" tmp_empirical/
        done <unique_genes_empirical.txt

        # Need to add -L to cp command when running on a compute nodes scratch space
        while read gene; do
          cp -rL "!{permuted}/phenotype=${gene}" tmp_permuted/
        done <unique_genes_permuted.txt

        run_susie_over_loci.R \
          --permuted tmp_permuted \
          --empirical tmp_empirical \
          --variant-reference !{variantReference} \
          --uncorrelated-genes unique_genes_permuted.txt \
          --bed-files !{bedFile.join(" ")}
        '''
}

process ExportResults {
  publishDir "${params.output}/finemapped", mode: 'move', overwrite: true

  input:
      path finemapped

  output:
      path "finemapped.results.tsv"

  shell:
      '''
      filter_finemapped_results.py \
      -i !{finemapped} \
      -s no_filter \
      -o finemapped.results.tsv
      '''
}
