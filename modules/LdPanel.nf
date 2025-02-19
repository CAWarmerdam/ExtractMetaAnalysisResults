#!/bin/bash nextflow


process SplitVariantSet {

    input:
        path variantReference
        val chunks

    output:
        path "chunks.txt"

    script:
        """
        export PYTHONUNBUFFERED="1"
        split_variant_set.py \
        --variant-reference ${variantReference} \
        --chunks ${chunks} \
        --output chunks.txt
        """
}


process GenerateLdPanel {
    publishDir "${params.output}", mode: 'copy', overwrite: true, pattern: 'chr=*/*.parquet'

    input:
        path parquetDataset
        path pcaFolder
        path variantReference
        path ldGenes
        tuple val(chromosome), val(min_variant_index), val(max_variant_index), val(chunk_nr)

    output:
        path "chr=${chromosome}/*parquet"

    script:
        """
        export PYTHONUNBUFFERED="1"
        prepare_ld_panel.py \
        --dataset-folder ${parquetDataset} \
        --genes-file ${ldGenes} \
        --variant-reference ${variantReference} \
        --pca-prefix "${pcaFolder}/coexp_pca_10496_genes_pruned" \
        --chromosome ${chromosome} \
        --min-max ${min_variant_index} ${max_variant_index} \
        --output-folder '.'
        """
}
