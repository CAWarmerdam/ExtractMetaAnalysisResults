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

process RunSusieFineMapping {
    //scratch true // Needs to be set to true!
    scratch '$TMPDIR'
    publishDir "${params.output}/finemapped", mode: 'copy', overwrite: true

    input:
        path empirical, stageAs: 'empirical'
        path permuted, stageAs: 'permuted'
        path variantReference
        path uncorrelatedGenes
        tuple val(chromosome), path(bedFile)
        val ld_type
        val max_i2

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
        cp -rL "!{permuted}/ld_panel_chr!{chromosome}.parquet" tmp_permuted/

        run_susie_over_loci.R \
          --ld tmp_permuted/ld_panel_chr!{chromosome}.parquet \
          --empirical tmp_empirical \
          --variant-reference !{variantReference} \
          --uncorrelated-genes unique_genes_permuted.txt \
          --bed-files !{bedFile.join(" ")} \
          --ld-type !{ld_type} \
          --max-i2 !{max_i2}

        rm -r tmp_empirical/
        rm -r tmp_permuted/
        '''
}

process RunCarmaFineMapping {
    //scratch true // Needs to be set to true!
    scratch '$TMPDIR'
    publishDir "${params.output}/finemapped", mode: 'copy', overwrite: true
    container null

    beforeScript "ml R/4.4.1-gfbf-2023b; ml GSL/2.7-GCC-13.2.0"

    input:
	path empirical, stageAs: 'empirical'
        path permuted, stageAs: 'permuted'
        path variantReference
        path uncorrelatedGenes
        tuple val(chromosome), path(bedFile)
        val ld_type
        val max_i2

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
        cp -rL "!{permuted}/ld_panel_chr!{chromosome}.parquet" tmp_permuted/

        run_carma_over_loci.R \
          --ld tmp_permuted/ld_panel_chr!{chromosome}.parquet \
          --empirical tmp_empirical \
          --variant-reference !{variantReference} \
          --uncorrelated-genes unique_genes_permuted.txt \
          --bed-files !{bedFile.join(" ")} \
          --ld-type !{ld_type} \
          --max-i2 !{max_i2}

        rm -r tmp_empirical/
        rm -r tmp_permuted/
        '''
}


process RunRSparseProFineMapping {
    //scratch true // Needs to be set to true!
    scratch '$TMPDIR'
    publishDir "${params.output}/finemapped", mode: 'copy', overwrite: true

    input:
    path empirical, stageAs: 'empirical'
        path permuted, stageAs: 'permuted'
        path variantReference
        path uncorrelatedGenes
        tuple val(chromosome), path(bedFile)
        val ld_type
        val max_i2

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
        cp -rL "!{permuted}/ld_panel_chr!{chromosome}.parquet" tmp_permuted/

        run_rsparsepro_over_loci.py \
          --ld tmp_permuted/ld_panel_chr!{chromosome}.parquet \
          --empirical tmp_empirical \
          --variant-reference !{variantReference} \
          --uncorrelated-genes unique_genes_permuted.txt \
          --bed-files !{bedFile.join(" ")} \
          --ld-type !{ld_type} \
          --max-i2 !{max_i2}

        rm -r tmp_empirical/
        rm -r tmp_permuted/
        '''
}


process RunCarmaSusieFineMapping {
    //scratch true // Needs to be set to true!
    scratch '$TMPDIR'
    publishDir "${params.output}/finemapped", mode: 'copy', overwrite: true

    input:
    path empirical, stageAs: 'empirical'
        path permuted, stageAs: 'permuted'
        path variantReference
        path uncorrelatedGenes
        tuple val(chromosome), path(bedFile)
        val ld_type
        val max_i2

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
        cp -rL "!{permuted}/ld_panel_chr!{chromosome}.parquet" tmp_permuted/

        run_carma_light_over_loci.py \
          --ld tmp_permuted/ld_panel_chr!{chromosome}.parquet \
          --empirical tmp_empirical \
          --variant-reference !{variantReference} \
          --uncorrelated-genes unique_genes_permuted.txt \
          --bed-files !{bedFile.join(" ")} \
          --ld-type !{ld_type} \
          --max-i2 !{max_i2}

        rm -r tmp_empirical/
        rm -r tmp_permuted/
        '''
}




process ExportResults {
  publishDir "${params.output}/finemapped", mode: 'move', overwrite: true

  input:
      path finemapped, name: 'finemapped.result.*.tsv'

  output:
      path "finemapped.results.tsv"

  shell:
      '''
      filter_finemapped_results.py \
      -i !{finemapped.join(" ")} \
      -s no_filter \
      -o finemapped.results.tsv
      '''
}
