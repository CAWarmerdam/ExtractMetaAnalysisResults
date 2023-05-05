#!/bin/bash nextflow


process SampleOverlapMatrix {

    publishDir "${params.output}/sample_overlap", mode: 'copy', overwrite: true

    input:
        path inclusionDir

    output:
        path "sample_overlap_matrix.csv"

    script:
        """
        sample_overlap_table.py ${inclusionOutput}/filter_logs_full.log ${inclusionOutput}/filter_logs_full.log
        """
}