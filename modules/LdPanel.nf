#!/bin/bash nextflow


process SplitVariantSet {

    input:
        path variantSet
        path variantReference
        val chunks

    output:
        path "chunks.txt"

    script:
        """
        export PYTHONUNBUFFERED="1"
        ld_panel/split_variant_set.py \
        --variant-set ${variantSet} \
        --variant-reference ${variantReference} \
        --chunks ${chunks} \
        --output chunks.txt
        """
}


process GenerateLdPanel {
    publishDir "${params.output}", mode: 'move', overwrite: true

    input:
        path parquetDataset
        path variantReference
        path ldGenes
        tuple val(chromosome), val(min_variant_index), val(max_variant_index), val(chunk_nr)

    output:
        path "chr*"

    script:
        """
        export PYTHONUNBUFFERED="1"
        ld_panel/prepare_ld_panel.py \
        --dataset-folder ${parquetDataset} \
        --genes-file ${ldGenes} \
        --variant-reference ${variantReference} \
        --chromosome ${chromosome} \
        --min-max ${min_variant_index} ${max_variant_index} \
        --output-folder '.'
        """
}
